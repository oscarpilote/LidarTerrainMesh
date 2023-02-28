#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "meshoptimizer/src/meshoptimizer.h"

#include "array.h"
#include "hash_table.h"
#include "math_utils.h"
#include "sys_utils.h"

#include "mesh.h"
#include "mesh_inria.h"
#include "mesh_ply.h"
#include "mesh_utils.h"
#include "vertex_table.h"

#include "las_normal.h"
#include "las_read.h"
#include "las_source.h"
#include "las_tree.h"

/******************************************************************************
 *
 * I. Args and filenames.
 *
 ******************************************************************************/

struct Cfg {
	int x0;
	int y0;
	int depth;
	int weight;
	char base_dir[128];
	char out_dir[128];
	int x0_base;
	int y0_base;
	bool verbose;
	bool simplify;
	bool optimize;
	bool encode;
	float hausd;
	float hgrad;
};

static int process_args(int argc, const char **argv, struct Cfg &cfg)
{
	if (argc < 3) {
		printf("Syntax : swiss_lidar x y [depth=10] [weight=8] "
		       "[base_dir=.]\n");
		return (-1);
	}
	cfg.x0 = atoi(argv[1]);
	cfg.y0 = atoi(argv[2]);
	cfg.depth = (argc >= 4) ? atoi(argv[3]) : 10;
	cfg.weight = (argc >= 5) ? atoi(argv[4]) : 8;
	if (argc >= 6) {
		if (strlen(argv[5]) > 127) {
			printf("Base_dir size overflow (max 127)\n");
			return (-1);
		}
		strncpy(cfg.base_dir, argv[5], 127);
	}
	if (argc >= 7) {
		if (strlen(argv[6]) > 127) {
			printf("Out_dir size overflow (max 127)\n");
			return (-1);
		}
		strncpy(cfg.out_dir, argv[6], 127);
	}
	cfg.x0_base = (argc >= 8) ? atoi(argv[7]) : 0;
	cfg.y0_base = (argc >= 9) ? atoi(argv[8]) : 0;
	cfg.verbose = (argc >= 10) ? atoi(argv[9]) : 1;
	cfg.simplify = (argc >= 11) ? atoi(argv[10]) : 0;
	cfg.optimize = (argc >= 12) ? atoi(argv[11]) : 0;
	cfg.encode = (argc >= 13) ? atoi(argv[12]) : 0;
	cfg.hausd = (argc >= 14) ? atof(argv[13]) : 0.15f;
	cfg.hgrad = (argc >= 15) ? atof(argv[14]) : 5.f;

	return (0);
}

static void print_cfg(struct Cfg &cfg)
{
	printf("\n");
	printf("Configuration :\n");
	printf("---------------\n");
	printf("Tile coords : %d %d\n", cfg.x0, cfg.y0);
	printf("Base coords : %d %d\n", cfg.x0_base, cfg.y0_base);
	printf("Data  dir   : %s\n", cfg.base_dir);
	printf("Output dir  : %s\n", cfg.out_dir);
	printf("Verbosity   : %d\n", cfg.verbose ? 1 : 0);
	if (cfg.simplify) {
		printf("Simplify    : hausd=%f , hgrad=%f\n", cfg.hausd,
		       cfg.hgrad);
	}
}

static char *las_filename(int x, int y, const char *base_dir)
{
	int base_len = 0;
	bool trailing_slash = true;
	if (base_dir) {
		base_len = strlen(base_dir);
		if (base_len && base_dir[base_len - 1] != '/') {
			trailing_slash = false;
			base_len += 1;
		}
	}
	char *fname = (char *)calloc((base_len + 16), sizeof(*fname));
	if (base_dir) {
		strcpy(fname, base_dir);
	}
	if (!trailing_slash) {
		fname[base_len - 1] = '/';
	}
	snprintf(fname + base_len, 16, "%4d_%4d.las", x, y);

	return fname;
}

static char *oriented_points_filename(int x, int y, const char *out_dir)
{
	int base_len = 0;
	bool trailing_slash = true;
	base_len = strlen(out_dir);
	if (base_len && out_dir[base_len - 1] != '/') {
		trailing_slash = false;
		base_len += 1;
	}
	char *fname = (char *)calloc((base_len + 22), sizeof(*fname));
	strcpy(fname, out_dir);
	if (!trailing_slash) {
		fname[base_len - 1] = '/';
	}
	snprintf(fname + base_len, 22, "%4d_%4d-points.ply", x, y);

	return fname;
}

static char *mesh_filename(int x, int y, const char *base_dir, const char *ext)
{
	int base_len = 0;
	bool trailing_slash = true;
	base_len = strlen(base_dir);
	if (base_len && base_dir[base_len - 1] != '/') {
		trailing_slash = false;
		base_len += 1;
	}
	char *fname = (char *)calloc((base_len + 32), sizeof(*fname));
	strcpy(fname, base_dir);
	if (!trailing_slash) {
		fname[base_len - 1] = '/';
	}
	snprintf(fname + base_len, 32, "%4d_%4d.%s", x, y, ext);

	return fname;
}

/******************************************************************************
 *
 * II. Functions related to creation of the oriented point set.
 *
 ******************************************************************************/

struct Transform {
	float scale;
	Vec3 shift;
};

static inline int get_source_idx(const TArray<int> &sources, int source_id)
{
	int n = sources.size;
	while (n--) {
		if (sources[n] == source_id)
			return (n);
	}
	return (-1);
}

static inline void rebase_point(struct LasPoint &p, int dx, int dy)
{
	p.x += dx * 100000;
	p.y += dy * 100000;
}

static inline bool filter_point(const struct LasPoint &p)
{
	/* 100m boundary buffer size */
	int bd = 10000;
	bool ret = (p.x > -bd) && (p.x < 100000 + bd) && (p.y > -bd) &&
		   (p.y < 100000 + bd) && (p.classification == 2);
	return (ret);
}

static int get_ids_and_counts(TArray<int> &source_ids,
			      TArray<int> &source_counts, int x0, int y0,
			      const char *base_dir)
{
	int file_count = 0;
	for (int i = 0; i < 9; ++i) {
		int dx = (i % 3) - 1;
		int dy = (i / 3) - 1;
		int x = x0 + dx;
		int y = y0 + dy;
		char *fname = las_filename(x, y, base_dir);
		struct LasFileInfo info;
		if (las_read_info(fname, info)) {
			free(fname);
			continue;
		}
		file_count += 1;

		char *raw_data = las_load_raw_data(fname, info, NULL);
		assert(raw_data);

		struct LasPoint p;
		for (size_t i = 0; i < info.point_num; ++i) {
			const char *raw = raw_data + i * info.point_size;
			p = las_read_point(raw, info.point_format);
			rebase_point(p, dx, dy);
			if (!filter_point(p))
				continue;
			int idx = get_source_idx(source_ids, p.source_id);
			if (idx >= 0) {
				source_counts[idx]++;
			} else {
				source_ids.push_back(p.source_id);
				source_counts.push_back(1);
			}
		}
		free(raw_data);
		free(fname);
	}
	return file_count;
}

static int accumulate_counts(TArray<int> &offsets, const TArray<int> &counts)
{
	assert(offsets.size == counts.size);
	offsets[0] = 0;
	for (size_t i = 1; i < offsets.size; ++i) {
		offsets[i] = offsets[i - 1] + counts[i - 1];
	}
	return offsets[offsets.size - 1] + counts[counts.size - 1];
}

static void fill_points(TArray<struct LasPoint> &points, const TArray<int> &ids,
			const TArray<int> &counts, TArray<int> &offsets, int x0,
			int y0, const char *base_dir)
{
	for (int i = 0; i < 9; ++i) {
		int dx = (i % 3) - 1;
		int dy = (i / 3) - 1;
		int x = x0 + dx;
		int y = y0 + dy;
		char *fname = las_filename(x, y, base_dir);
		struct LasFileInfo info;
		if (las_read_info(fname, info)) {
			free(fname);
			continue;
		}

		char *raw_data = las_load_raw_data(fname, info, NULL);
		assert(raw_data);

		struct LasPoint p;
		for (size_t i = 0; i < info.point_num; ++i) {
			const char *raw = raw_data + i * info.point_size;
			p = las_read_point(raw, info.point_format);
			rebase_point(p, dx, dy);
			if (!filter_point(p))
				continue;
			int idx = get_source_idx(ids, p.source_id);
			assert(idx >= 0);
			p.source_idx = idx;
			points[offsets[idx]] = p;
			/* We will settle them back before exit */
			offsets[idx]++;
		}
		free(raw_data);
		free(fname);
	}
	/* Reset offsets to their correct values */
	for (size_t i = 0; i < offsets.size; ++i) {
		offsets[i] -= counts[i];
	}
}

static int read_and_filter_data(size_t &source_num,
				TArray<struct LasPoint> &points,
				const struct Cfg &cfg)
{
	TArray<int> ids;
	TArray<int> counts;
	get_ids_and_counts(ids, counts, cfg.x0, cfg.y0, cfg.base_dir);
	source_num = ids.size;
	if (!source_num)
		return (-1);

	TArray<int> offsets(source_num);
	int point_num = accumulate_counts(offsets, counts);
	points.clear();
	points.resize(point_num);
	fill_points(points, ids, counts, offsets, cfg.x0, cfg.y0, cfg.base_dir);

	return (0);
}

static int send_points_to_unit_cube(const TArray<struct LasPoint> &points,
				    Vec3 *pos, struct Transform &t)
{
	int min_z = INT_MAX, max_z = -INT_MAX;
	for (size_t i = 0; i < points.size; ++i) {
		min_z = MIN(min_z, points[i].z);
		max_z = MAX(max_z, points[i].z);
	}
	float span = (max_z - min_z) * 0.01f;
	printf("(Altitude span : %.0f meters)\n", span);
	float n;
	if (span < 1083.33f) {
		n = 6;
	} else if (span < 1625.f) {
		n = 2;
	} else {
		return (-1);
	}
	t.scale = n / (n + 2) * 1e-5;
	t.shift.x = t.shift.y = 1.f / (n + 2);
	float mean = 0.5 * (min_z + max_z) * t.scale;
	/* Round shift.z so that octree boxes match vertically to neighboors */
	t.shift.z = 0.5 - round(16 * mean) * 0.0625;
	assert(min_z * t.scale + t.shift.z >= 0.0625);
	assert(max_z * t.scale + t.shift.z <= 1 - 0.0625);
	for (size_t i = 0; i < points.size; ++i) {
		pos[i].x = points[i].x * t.scale + t.shift.x;
		pos[i].y = points[i].y * t.scale + t.shift.y;
		pos[i].z = points[i].z * t.scale + t.shift.z;
	}
	return (0);
}

static void nml_cb(float progress)
{
	printf("\rEval. normal directions    : %.0f %%", progress);
	fflush(stdout);
}

static int build_oriented_point_set(Mesh &mesh, MBuf &data,
				    struct Transform &transf,
				    const struct Cfg &cfg)
{
	/* Read data */
	size_t source_num;
	TArray<struct LasPoint> points;
	if (read_and_filter_data(source_num, points, cfg)) {
		return (-1);
	}
	size_t point_num = points.size;
	printf("Total Lidar sources        : %zu\n", source_num);
	printf("Total filtered las points  : %zu\n", point_num);

	/* Derive fligh lines azimut from data. Used later for normals */
	printf("Reconstructing flightlines :");
	TArray<struct SourceStat> stats(source_num);
	las_stat_sources(points, stats);
	TArray<struct SourceFlightLine> fls(source_num);
	double scale[3] = {1., 1., 1.};
	int valid = las_approx_flight_lines(points, scale, stats, fls);
	printf(" %d valid\n", valid);

	printf("Set positions & transform  : ");
	/* Rescale and offset positions into buffer */
	mesh.vertex_count = point_num;
	data.vtx_attr = VtxAttr::PN;
	data.reserve_vertices(point_num);
	if (send_points_to_unit_cube(points, data.positions, transf)) {
		printf("Altitude span too large !\n");
		return (-1);
	}

	/* Build kd-tree for filtered point cloud. */
	struct LasCoords coords = {data.positions, point_num};
	LasTree tree(3, coords, 10);

	/* Build normals using geometry and scanlines*/
	printf("Eval. normal directions    :  ");
	TArray<float> qual(point_num);
	estim_unoriented_nml(data.positions, point_num, data.normals, qual.data,
			     tree, nml_cb);
	printf("\n");

	printf("Eval. normal orientations  :\n");
	TArray<EOrient> oriented(point_num, ENone);

	size_t unset;

	// TODO: Must take quality into account for z pass
	// unset = orient_normals_with_z(nml, oriented);
	// printf("    Oriented after positive z pass : %.1f %%\n",
	//        (nml.size - unset) * 100.f / nml.size);

	unset =
	    orient_nml_with_scan(points.data, point_num, fls.data, source_num,
				 qual.data, data.normals, oriented.data);
	printf("    Oriented after scanlines pass : %.1f %%\n",
	       (point_num - unset) * 100.f / point_num);

	float progress = 1.f;
	int pass = 1;
	while (unset && progress > 0.01f) {
		size_t newly_set =
		    propagate_nml_once(data.positions, point_num, tree,
				       qual.data, data.normals, oriented.data);
		progress = (float)newly_set / unset;
		unset -= newly_set;
		printf("    Oriented after propagate pass nbr %d : %.1f %%\n",
		       pass++, (point_num - unset) * 100.f / point_num);
	}

	/* Skip points with no orientation found */
	{
		size_t new_num = 0;
		for (size_t i = 0; i < point_num; ++i) {
			if (oriented[i] != ENone) {
				data.positions[new_num] = data.positions[i];
				data.normals[new_num] = data.normals[i];
				new_num++;
			}
		}
		printf("    %.1f %% vertices discarded due to lack of clear "
		       "orientation.\n",
		       unset * 100.f / point_num);
		mesh.vertex_count = new_num;
	}

	/* Add dummy points for bounding box to perfectly match unit cube */
	data.reserve_vertices(mesh.vertex_count + 2);
	data.positions[mesh.vertex_count] = Vec3{0.f, 0.f, 0.f};
	data.normals[mesh.vertex_count++] = Vec3{0.f, 0.f, 1.f};
	data.positions[mesh.vertex_count] = Vec3{1.f, 1.f, 1.f};
	data.normals[mesh.vertex_count++] = Vec3{0.f, 0.f, 1.f};

	return (0);
}

/******************************************************************************
 *
 * III. Functions related to surface mesh (post)processing.
 *
 ******************************************************************************/

static int launch_poisson_recon(Mesh &mesh, MBuf &data, const struct Cfg &cfg)
{
	char *recon_in = oriented_points_filename(cfg.x0, cfg.y0, cfg.out_dir);

	/* TODO: Avoid going to disk */
	write_ply(recon_in, mesh, data);

	char *recon_out = mesh_filename(cfg.x0, cfg.y0, cfg.out_dir, "ply");
	const char *format =
	    "poissonrecon --in %s --out %s --scale 1.0 --depth %d "
	    "--pointWeight %d --threads 8 %s";
	unsigned short len =
	    strlen(recon_in) + strlen(recon_out) + strlen(format) + 16;

	char *cmd = (char *)calloc(len, sizeof(*cmd));
	snprintf(cmd, len, format, recon_in, recon_out, cfg.depth, cfg.weight,
		 cfg.verbose ? "--verbose" : "");
	int ret = system(cmd);

	snprintf(cmd, len, "rm %s", recon_in);
	system(cmd);

	if (!ret) {
		load_ply(mesh, data, recon_out);
	}

	free(cmd);
	free(recon_out);
	free(recon_in);

	return (ret);
}

static int postprocess_surface_mesh(Mesh &mesh, MBuf &data,
				    const struct Cfg &cfg,
				    const struct Transform &transf)
{
	/* TODO Avoid this copy and swap of mesh/buffers, this requires
	 * implementing vertex swaps in join_mesh_from_indices */
	Mesh mesh2;
	MBuf data2;
	data2.vtx_attr = VtxAttr::POS;
	data2.reserve_indices(mesh.index_count);
	data2.reserve_vertices(mesh.vertex_count);

	/* Remove triangles in buffered part */
	for (size_t i = 0; i < mesh.index_count / 3; ++i) {
		uint32_t i0 = data.indices[3 * i + 0];
		uint32_t i1 = data.indices[3 * i + 1];
		uint32_t i2 = data.indices[3 * i + 2];
		Vec3 p0 = data.positions[i0];
		Vec3 p1 = data.positions[i1];
		Vec3 p2 = data.positions[i2];
		Vec3 bary = (p0 + p1 + p2) * (1 / 3.f);
		if (bary.x < transf.shift.x || bary.x > 1.f - transf.shift.x ||
		    bary.y < transf.shift.x || bary.y > 1.f - transf.shift.y) {
			continue;
		}
		data2.indices[mesh2.index_count++] = i0;
		data2.indices[mesh2.index_count++] = i1;
		data2.indices[mesh2.index_count++] = i2;
	}

	printf("A total of %.1f MTri after buffer cut.\n",
	       1e-6 * mesh2.index_count / 3);

	/* Inverse transform points + shift relative to base */
	/* TODO : should be eventually removed and use scene
	 *        object placement and scaling instead
	 */
	Vec3 base_shift{(cfg.x0 - cfg.x0_base) * 1000.f,
			(cfg.y0 - cfg.y0_base) * 1000.f, 0};
	float scale = (1 / (100 * transf.scale)); /* in meters */
	for (size_t i = 0; i < mesh.vertex_count; ++i) {
		data2.positions[i] = data.positions[i] - transf.shift;
		data2.positions[i] *= scale;
		data2.positions[i] += base_shift;
	}
	mesh2.vertex_count = mesh.vertex_count;

	/* Compact mesh (useless vertices and degenerate triangles removed) */
	mesh.clear();
	data.update_vtx_attr(VtxAttr::POS);
	VertexTable vtx_table(mesh2.vertex_count, &data, data.vtx_attr);
	join_mesh_from_indices(mesh, data, mesh2, data2, vtx_table, NULL);
	skip_degenerate_tris(mesh, data);

	printf("A total of %.1f MTri after compacting.\n",
	       1e-6 * mesh.index_count / 3);

	return (0);
}

static int launch_mesh_simplification(Mesh &mesh, MBuf &data,
				      const struct Cfg &cfg)
{
	char *fin = mesh_filename(cfg.x0, cfg.y0, cfg.out_dir, "mesh");
	char *fout = mesh_filename(cfg.x0, cfg.y0, cfg.out_dir, "out.mesh");
	char *fsol = mesh_filename(cfg.x0, cfg.y0, cfg.out_dir, "out.sol");

	/* TODO Avoid going to disk */
	write_inria(fin, mesh, data);

	const char *format =
	    "mmgs -in %s -out %s -hausd %.1f -hgrad %.1f -nr -v %d";
	int len = strlen(fin) + strlen(fout) + strlen(format) + 8;

	char *cmd = (char *)calloc(len, sizeof(*cmd));
	snprintf(cmd, len, format, fin, fout, cfg.hausd, cfg.hgrad,
		 cfg.verbose ? 1 : -1);
	int ret = system(cmd);

	snprintf(cmd, len, "rm %s", fin);
	system(cmd);

	if (!ret) {
		read_inria(mesh, data, fout);
		snprintf(cmd, len, "rm %s", fout);
		system(cmd);
		snprintf(cmd, len, "rm %s", fsol);
		system(cmd);
	}

	free(cmd);
	free(fsol);
	free(fout);
	free(fin);

	return (ret);
}

int optimize_mesh(Mesh &mesh, MBuf &data, const Cfg &cfg)
{
	uint32_t index_count = mesh.index_count;
	uint32_t *indices = data.indices + mesh.index_offset;
	uint32_t vertex_count = mesh.vertex_count;
	float *vertices = (float *)(data.positions + mesh.vertex_offset);
	size_t vertex_size = sizeof(Vec3);
	meshopt_optimizeVertexCache(indices, indices, mesh.index_count, mesh.vertex_count);
	/* TODO This only work for POS only meshes, use
	 * meshopt_optimizeVertexFetchRemap instead
	 * */
	meshopt_optimizeVertexFetch(vertices, indices, index_count, vertices, vertex_count, vertex_size);
	return (0);
}

int quantize_encode_mesh(Mesh &mesh, MBuf &data, const Cfg &cfg, const Transform &transf)
{
	uint32_t index_count = mesh.index_count;
	uint32_t vertex_count = mesh.vertex_count;
	Vec3 base_shift{(cfg.x0 - cfg.x0_base) * 1000.f,
			(cfg.y0 - cfg.y0_base) * 1000.f, 0};
	float invscale = (100 * transf.scale);
	TArray<TVec3<uint16_t>> qpos(vertex_count);
	for (size_t i = 0; i < vertex_count; ++i) {
		Vec3 cubepos = data.positions[i];
		cubepos -= base_shift;
		cubepos *= invscale;
		cubepos += transf.shift;
		qpos[i].x = cubepos.x * ((1 << 16) - 1);
		qpos[i].y = cubepos.y * ((1 << 16) - 1);
		qpos[i].z = cubepos.z * ((1 << 16) - 1);
	}
	uint32_t *indices = data.indices + mesh.index_offset;
	//void *vertices = (float *)(data.positions + mesh.vertex_offset);
	//size_t vertex_size = sizeof(Vec3);
	void *vertices = qpos.data;
	size_t vertex_size = sizeof(TVec3<uint16_t>);
	TArray<uint8_t> vbuf(meshopt_encodeVertexBufferBound(vertex_count, vertex_size));
	vbuf.resize(meshopt_encodeVertexBuffer(&vbuf[0], vbuf.size, vertices, vertex_count, vertex_size));
	TArray<uint8_t> ibuf(meshopt_encodeIndexBufferBound(index_count, vertex_count));
	ibuf.resize(meshopt_encodeIndexBuffer(&ibuf[0], ibuf.size, indices, index_count));

	char *fname = mesh_filename(cfg.x0, cfg.y0, cfg.out_dir, "meshopt");
	FILE *f = fopen(fname, "wb");
	int ret = (f == NULL)
		|| fwrite(vbuf.data, vbuf.size, 1, f) != 1 
		|| fwrite(ibuf.data, ibuf.size, 1, f) != 1
		? -1 : 0;
	fclose(f);
	free(fname);
	return (ret);
}


int write_final_mesh(const Mesh &mesh, const MBuf &data, const struct Cfg &cfg)
{

	const char *ext = cfg.simplify ? "mmgs.ply" : "ply";

	char *fname = mesh_filename(cfg.x0, cfg.y0, cfg.out_dir, ext);
	write_ply(fname, mesh, data);
	free(fname);

	return (0);
}

/******************************************************************************
 *
 * IV. Main.
 *
 ******************************************************************************/

int main(int argc, char **argv)
{
	struct Cfg cfg;
	Mesh mesh;
	MBuf data;

	/* Process command line arguments */
	if (process_args(argc, (const char **)argv, cfg)) {
		return (-1);
	}
	printf("\n------ Dealing with %04d %04d ------\n", cfg.x0, cfg.y0);
	print_cfg(cfg);

	/* Gather neighboring las files to build a 100m buffer of the
	 * target tile, and then compute combined position + normal point
	 * set from it. Normal orientation reconstruction combines geometric
	 * and lidar source statistics which are computed on their own.
	 */
	printf("\n");
	printf("I. Building oriented point set :\n");
	printf("--------------------------------\n");
	struct Transform transf = {1, Vec3::Zero};
	if (build_oriented_point_set(mesh, data, transf, cfg)) {
		printf("No data found.\n");
		return (-1);
	}

	printf("\n");
	printf("II. Building surface mesh from point set :\n");
	printf("------------------------------------------\n");
	/* Launch 3rd party screened Poisson Reconstruction.
	 */
	if (launch_poisson_recon(mesh, data, cfg)) {
		printf("Error in Poisson reconstruction\n");
		return (-1);
	}
	printf("Poisson reconstruction completed : a total of %.1f MTri and "
	       "%.1f MVert.\n",
	       1e-6 * mesh.index_count / 3, 1e-6 * mesh.vertex_count);

	printf("\n");
	printf("III. Postprocessing surface mesh.\n");
	printf("---------------------------------\n");
	/* Remove buffered boundary and apply some first cleaning */
	if (postprocess_surface_mesh(mesh, data, cfg, transf)) {
		printf("Error in Postprocessing\n");
		return (-1);
	}

	/* Apply simplification filter */
	if (cfg.simplify && launch_mesh_simplification(mesh, data, cfg)) {
		printf("Error in Simplification\n");
		return (-1);
	}
	/* Apply optimization filter */
	if (cfg.optimize && optimize_mesh(mesh, data, cfg)) {
		printf("Error in Mesh Optimization\n");
		return (-1);
	}

	/* Save final mesh */
	if (cfg.encode && quantize_encode_mesh(mesh, data, cfg, transf)) {
		printf("Error in mesh encoding\n");
		return (-1);
	}
	
	write_final_mesh(mesh, data, cfg);

	printf("\n------ Finished with %04d %04d ------\n", cfg.x0, cfg.y0);
	return (0);
}

/******************************************************************************/
