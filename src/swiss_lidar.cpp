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
#include "mesh_remesh.h"
#include "mesh_ply.h"
#include "mesh_utils.h"
#include "vertex_table.h"

#include "las_normal.h"
#include "las_read.h"
#include "las_source.h"

#include "swiss_raster.h"

#define DEBUG_NML_CONFIDENCE 1

/******************************************************************************
 *
 * I. Args and filenames.
 *
 ******************************************************************************/

struct Cfg {
	int x0;
	int y0;
	char base_dir[128];
	char out_dir[128];
	int depth;
	int weight;
	float simp_error;
	float hausd;
	float hgrad;
	bool ani;
	bool clean;
	bool verbose;
	bool optimize;
	bool encode;
	bool nml_confidence;
	int x0_base;
	int y0_base;
};

static int process_args(int argc, const char **argv, struct Cfg &cfg)
{
	if (argc < 3) {
		printf("Syntax : swiss_lidar x y [...]\n");
		return (-1);
	}
	cfg.x0 = atoi(argv[1]);
	cfg.y0 = atoi(argv[2]);
	if (argc >= 4) {
		if (strlen(argv[3]) > 127) {
			printf("Base_dir size overflow (max 127)\n");
			return (-1);
		}
		strncpy(cfg.base_dir, argv[3], 127);
	} else {
		strcpy(cfg.base_dir, ".");
	}
	if (argc >= 5) {
		if (strlen(argv[4]) > 127) {
			printf("Out_dir size overflow (max 127)\n");
			return (-1);
		}
		strncpy(cfg.out_dir, argv[4], 127);
	} else {
		strcpy(cfg.base_dir, ".");
	}
	cfg.depth = (argc >= 6) ? atoi(argv[5]) : 10;
	cfg.weight = (argc >= 7) ? atoi(argv[6]) : 8;
	cfg.simp_error = (argc >= 8) ? atof(argv[7]) : 1e-4;
	cfg.hausd = (argc >= 9) ? atof(argv[8]) : 0.15f;
	cfg.hgrad = (argc >= 10) ? atof(argv[9]) : 5.f;
	cfg.ani = (argc >= 11) ? atoi(argv[10]) : 0;
	cfg.clean = (argc >= 12) ? atoi(argv[11]) : 0;
	cfg.verbose = (argc >= 13) ? atoi(argv[12]) : 0;
	cfg.optimize = (argc >= 14) ? atoi(argv[13]) : 1;
	cfg.encode = (argc >= 15) ? atoi(argv[14]) : 0;
	cfg.nml_confidence = (argc >= 16) ? atoi(argv[15]) : 1;
	cfg.x0_base = (argc >= 17) ? atoi(argv[16]) : cfg.x0;
	cfg.y0_base = (argc >= 18) ? atoi(argv[17]) : cfg.y0;

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
	printf("Simplify    : meshopt=%f,  hausd=%f , hgrad=%f\n",
	       cfg.simp_error, cfg.hausd, cfg.hgrad);
}

static char *get_filename(int x, int y, const char *dir, const char *ext)
{
	int base_len = 0;
	bool trailing_slash = true;
	base_len = strlen(dir);
	if (base_len && dir[base_len - 1] != '/') {
		trailing_slash = false;
		base_len += 1;
	}
	char *fname = (char *)calloc((base_len + 32), sizeof(*fname));
	strcpy(fname, dir);
	if (!trailing_slash) {
		fname[base_len - 1] = '/';
	}
	snprintf(fname + base_len, 32, "%4d_%4d.%s", x, y, ext);

	return fname;
}

/******************************************************************************
 *
 * II. Utility fonctions.
 *
 ******************************************************************************/

static int write_mesh(const Mesh &mesh, const MBuf &data, const struct Cfg &cfg,
		      const char *ext)
{

	char *fname = get_filename(cfg.x0, cfg.y0, cfg.out_dir, ext);
	write_ply(fname, mesh, data);
	free(fname);

	return (0);
}

static void compact_mesh(Mesh &mesh, MBuf &data)
{
	/* TODO Avoid this copy and swap of mesh/buffers, this requires
	 * implementing vertex swaps in join_mesh_from_indices */
	Mesh mesh2 = mesh;
	MBuf data2;
	data2.vtx_attr = VtxAttr::POS;
	data2.reserve_indices(mesh.index_count);
	data2.reserve_vertices(mesh.vertex_count);

	copy_indices(data2, 0, data, 0, mesh.index_count);
	copy_vertices(data2, 0, data, 0, mesh.vertex_count);

	/* Compact mesh (useless vertices and degenerate triangles removed) */
	mesh.clear();
	data.update_vtx_attr(VtxAttr::POS);
	VertexTable vtx_table(mesh2.vertex_count, &data, data.vtx_attr);
	join_mesh_from_indices(mesh, data, mesh2, data2, vtx_table, NULL);
	skip_degenerate_tris(mesh, data);

	printf("A total of %d (%.2f M) Tri after compacting.\n",
	       mesh.index_count / 3, 1e-6 * mesh.index_count / 3);
}

#if 0
static void quantize_mesh(const Mesh &mesh, MBuf &data, float h)
{
	if (h <= 0)
		return;
	float *pos = (float *)(data.positions + mesh.vertex_offset);
	float invh = 1. / h;
	for (size_t i = 0; i < mesh.vertex_count; ++i) {
		pos[3 * i + 0] = h * roundf(pos[3 * i + 0] * invh);
		pos[3 * i + 1] = h * roundf(pos[3 * i + 1] * invh);
		pos[3 * i + 2] = h * roundf(pos[3 * i + 2] * invh);
	}
}
#endif

#if 0
static void clusterize_mesh(const Mesh &mesh, MBuf &data, float h)
{
	if (h <= 0)
		return;

	const float sqr_h = h * h;
	const uint32_t todo = ~0u;

	TArray<uint32_t> remap(mesh.vertex_count, todo);

	const float *pos = (const float *)(data.positions + mesh.vertex_offset);
	KdCoords<3> kdcoords{pos, mesh.vertex_count};
	KdTree<3> kdtree(3, kdcoords, 10);
	const size_t probes = 10;
	TArray<unsigned> knn_idx(probes);
	TArray<float> sqr_dist(probes);
	for (size_t i = 0; i < mesh.vertex_count; ++i) {
		if (remap[i] != todo)
			continue;
		remap[i] = i;
		const float *query = &pos[3 * i];
		int found =
		    kdtree.knnSearch(query, probes, &knn_idx[0], &sqr_dist[0]);
		for (int j = 0; j < found; ++j) {
			if (sqr_dist[j] <= sqr_h)
				remap[knn_idx[j]] = i;
		}
	}
	uint32_t *indices = data.indices + mesh.index_offset;
	for (size_t i = 0; i < mesh.index_count; ++i) {
		indices[i] = remap[indices[i]];
	}
}
#endif

static void optimize_mesh(Mesh &mesh, MBuf &data)
{
	uint32_t index_count = mesh.index_count;
	uint32_t *indices = data.indices + mesh.index_offset;
	uint32_t vertex_count = mesh.vertex_count;
	float *vertices = (float *)(data.positions + mesh.vertex_offset);
	size_t vertex_size = sizeof(Vec3);
	meshopt_optimizeVertexCache(indices, indices, mesh.index_count,
				    mesh.vertex_count);
	/* TODO This only work for POS only meshes, use
	 * meshopt_optimizeVertexFetchRemap instead
	 * */
	mesh.vertex_count =
	    meshopt_optimizeVertexFetch(vertices, indices, index_count,
					vertices, vertex_count, vertex_size);
}

struct Transform {
	float scale;
	Vec3 shift;
};

int write_transform(const struct Transform &t, const struct Cfg &cfg)
{
	char *fname = get_filename(cfg.x0, cfg.y0, cfg.out_dir, "transf");
	FILE *f = fopen(fname, "w");
	if (!f) {
		free(fname);
		return (-1);
	}
	fprintf(f, "Scale %g\n", t.scale);
	fprintf(f, "Offset %g %g %g\n", t.shift.x, t.shift.y, t.shift.z);
	fclose(f);
	free(fname);
	return (0);
}

int read_transform(struct Transform &t, const struct Cfg &cfg)
{
	char *fname = get_filename(cfg.x0, cfg.y0, cfg.out_dir, "transf");
	FILE *f = fopen(fname, "r");
	if (!f) {
		free(fname);
		return (-1);
	}
	fscanf(f, "Scale %g\n", &t.scale);
	fscanf(f, "Offset %g %g %g\n", &t.shift.x, &t.shift.y, &t.shift.z);
	fclose(f);
	free(fname);
	return (0);
}

/******************************************************************************
 *
 * II. Functions related to creation of the oriented point set.
 *
 ******************************************************************************/

static inline int get_source_idx(const TArray<int> &sources, int source_id)
{
	int n = sources.size;
	while (n--) {
		if (sources[n] == source_id)
			return (n);
	}
	return (-1);
}

static inline void rebase_las_point(struct LasPoint &p, int dx, int dy)
{
	p.x += dx * 100000;
	p.y += dy * 100000;
}

static inline bool filter_las_point(const struct LasPoint &p)
{
	/* 100m boundary buffer size */
	int bd = 10000;
	bool ret = (p.x > -bd) && (p.x < 100000 + bd) && (p.y > -bd) &&
		   (p.y < 100000 + bd) && (p.classification == 2) ;
	return (ret);
}

static bool filter_raster_point(const struct LasPoint &p, const float *alt,
				int px, int py, int nx, int ny)
{
	/* 100m boundary buffer size */
	int bd = 10000;
	if ((p.x < -bd) || (p.x > 100000 + bd) || (p.y < -bd) ||
	    (p.y > 100000 + bd)) {
		return false;
	}

	/* Keep only points in non too steep terrain */
	int r = px + 1 < nx ? px + 1 : px;
	int l = px > 0 ? px - 1 : 0;
	int u = py + 1 < ny ? py + 1 : py;
	int d = py > 0 ? py - 1 : 0;

	float alt_p = alt[py * nx + px];
	float alt_r = alt[py * nx + r];
	float alt_l = alt[py * nx + l];
	float alt_u = alt[u * nx + px];
	float alt_d = alt[d * nx + px];

	float min_x, max_x, min_y, max_y;
	min_x = max_x = alt_p;
	min_x = MIN(min_x, alt_r);
	min_x = MIN(min_x, alt_l);
	max_x = MAX(max_x, alt_r);
	max_x = MAX(max_x, alt_l);
	min_y = max_y = alt_p;
	min_y = MIN(min_y, alt_u);
	min_y = MIN(min_y, alt_d);
	max_y = MAX(max_y, alt_u);
	max_y = MAX(max_y, alt_d);

	/* kind of... */
	float slope_x = (max_x - min_x) * nx / (2 * 1000);
	float slope_y = (max_y - min_y) * ny / (2 * 1000);
	float slope = MAX(slope_x, slope_y);

	return slope < 1;
}

static int get_las_counts(TArray<int> &source_ids, TArray<int> &source_counts,
			  int x0, int y0, const char *base_dir)
{
	int file_count = 0;
	for (int i = 0; i < 9; ++i) {
		int dx = (i % 3) - 1;
		int dy = (i / 3) - 1;
		int x = x0 + dx;
		int y = y0 + dy;
		char *fname = get_filename(x, y, base_dir, "las");
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
			rebase_las_point(p, dx, dy);
			if (!filter_las_point(p))
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

static void fill_las_points(TArray<struct LasPoint> &points,
			    const TArray<int> &ids, const TArray<int> &counts,
			    TArray<int> &offsets, int x0, int y0,
			    const char *base_dir)
{
	for (int i = 0; i < 9; ++i) {
		int dx = (i % 3) - 1;
		int dy = (i / 3) - 1;
		int x = x0 + dx;
		int y = y0 + dy;
		char *fname = get_filename(x, y, base_dir, "las");
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
			rebase_las_point(p, dx, dy);
			if (!filter_las_point(p))
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

static size_t read_and_filter_las_data(TArray<struct LasPoint> &points,
				       TArray<struct SourceFlightLine> &fls,
				       const struct Cfg &cfg,
				       bool verbose = true)
{
	TArray<int> ids;
	TArray<int> counts;
	get_las_counts(ids, counts, cfg.x0, cfg.y0, cfg.base_dir);
	size_t las_source_num = ids.size;
	if (!las_source_num) {
		if (verbose)
			printf("No LIDAR data found.\n");
		return (0);
	}
	TArray<int> offsets(las_source_num);
	int point_num = accumulate_counts(offsets, counts);
	points.clear();
	points.resize(point_num);
	fill_las_points(points, ids, counts, offsets, cfg.x0, cfg.y0,
			cfg.base_dir);
	/* Derive fligh lines azimut from data. Used later for normals
	 */
	if (verbose)
		printf("Reconstructing flightlines :");
	TArray<struct SourceStat> stats(las_source_num);
	las_stat_sources(points, stats);
	fls.resize(las_source_num);
	double scale[3] = {1., 1., 1.};
	int valid = las_approx_flight_lines(points, scale, stats, fls);
	if (verbose)
		printf(" found %d valid out of %zu sources\n", valid,
		       las_source_num);

	return points.size;
}

static size_t get_raster_point_count(int x0, int y0, const char *base_dir)
{
	size_t raster_point_count = 0;
	for (int i = 0; i < 9; ++i) {
		int dx = (i % 3) - 1;
		int dy = (i / 3) - 1;
		int x = x0 + dx;
		int y = y0 + dy;
		char *fname = get_filename(x, y, base_dir, "tif");
		int nx;
		int ny;
		if (read_swiss_geotiff(NULL, &nx, &ny, fname, true) || !nx ||
		    !ny) {
			free(fname);
			continue;
		}

		TArray<float> altitudes(nx * ny);
		read_swiss_geotiff(altitudes.data, &nx, &ny, fname);
		free(fname);
		for (int py = 0; py < ny; ++py) {
			for (int px = 0; px < nx; ++px) {
				struct LasPoint p;
				p.x = (2 * px + 1) * 50 * 1000 / nx;
				p.y = (2 * (ny - py) - 1) * 50 * 1000 / ny;
				p.z = (int)(100 * altitudes[nx * py + px]);
				rebase_las_point(p, dx, dy);
				if (filter_raster_point(p, altitudes.data, px,
							py, nx, ny))
					++raster_point_count;
			}
		}
	}
	return (raster_point_count);
}

static size_t fill_raster_points(TArray<struct LasPoint> &points, size_t offset,
				 int dummy_source_idx, int x0, int y0,
				 const char *base_dir)
{
	size_t point_idx = offset;
	for (int i = 0; i < 9; ++i) {
		int dx = (i % 3) - 1;
		int dy = (i / 3) - 1;
		int x = x0 + dx;
		int y = y0 + dy;
		char *fname = get_filename(x, y, base_dir, "tif");
		int nx;
		int ny;
		if (read_swiss_geotiff(NULL, &nx, &ny, fname, true) || !nx ||
		    !ny) {
			free(fname);
			continue;
		}

		TArray<float> altitudes(nx * ny);
		read_swiss_geotiff(altitudes.data, &nx, &ny, fname);
		free(fname);
		for (int py = 0; py < ny; ++py) {
			for (int px = 0; px < nx; ++px) {
				struct LasPoint p;
				p.x = (2 * px + 1) * 50 * 1000 / nx;
				p.y = (2 * (ny - py) - 1) * 50 * 1000 / ny;
				p.z = (int)(100 * altitudes[nx * py + px]);
				p.source_idx = dummy_source_idx;
				rebase_las_point(p, dx, dy);
				if (filter_raster_point(p, altitudes.data, px,
							py, nx, ny))
					points[point_idx++] = p;
			}
		}
	}
	return (point_idx - offset);
}

static size_t read_and_filter_raster_data(TArray<struct LasPoint> &points,
					  TArray<struct SourceFlightLine> &fls,
					  const struct Cfg &cfg)
{
	size_t raster_point_count =
	    get_raster_point_count(cfg.x0, cfg.y0, cfg.base_dir);

	if (!raster_point_count)
		return 0;

	size_t offset = points.size;
	points.reserve(points.size + raster_point_count);
	points.resize(points.size + raster_point_count);
	size_t dummy_source_idx = fls.size;
	fill_raster_points(points, offset, dummy_source_idx, cfg.x0, cfg.y0,
			   cfg.base_dir);
	fls.resize(fls.size + 1);
	fls[fls.size - 1].is_valid = false;

	return raster_point_count;
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
	t.scale = (float)n / (n + 2);
	t.shift.x = t.shift.y = 1.f / (n + 2);
	float mean = 0.5 * (min_z + max_z) * 1e-5 * t.scale;
	/* Round shift.z so that octree boxes match vertically to neighboors */
	t.shift.z = 0.5 - round(16 * mean) * 0.0625;
	assert(min_z * 1e-5 * t.scale + t.shift.z >= 0.0625);
	assert(max_z * 1e-5 * t.scale + t.shift.z <= 1 - 0.0625);
	for (size_t i = 0; i < points.size; ++i) {
		pos[i].x = points[i].x * 1e-5 * t.scale + t.shift.x;
		pos[i].y = points[i].y * 1e-5 * t.scale + t.shift.y;
		pos[i].z = points[i].z * 1e-5 * t.scale + t.shift.z;
	}
	return (0);
}

static void nml_cb(float progress)
{
	printf("\rEval. normal directions    : %.0f %%", progress);
	fflush(stdout);
}

/* Gather neighboring las files to build a 100m buffer of the
 * target tile, and then compute combined position + normal point
 * set from it. Normal orientation reconstruction combines geometric
 * and lidar source statistics which are computed on their own.
 */
static int build_oriented_point_set(const struct Cfg &cfg)
{

	/* Re-use existing output ? */
	{
		char *recon_in =
		    get_filename(cfg.x0, cfg.y0, cfg.out_dir, "points.ply");
		FILE *f;
		if ((f = fopen(recon_in, "rb")) != NULL) {
			printf("Using cached data in %s\n", recon_in);
			fclose(f);
			free(recon_in);
			return 0;
		}
	}

	/* Read data */
	TArray<struct LasPoint> points;
	TArray<struct SourceFlightLine> fls;
	// HACK
	// size_t las_num = 0;
	size_t las_num = read_and_filter_las_data(points, fls, cfg);
	size_t ras_num = read_and_filter_raster_data(points, fls, cfg);

	printf("Total Lidar points used    : %zu\n", las_num);
	printf("Total Raster points used   : %zu\n", ras_num);

	printf("Set positions & transform  : ");
	/* Rescale and offset positions into buffer */
	Mesh mesh;
	mesh.vertex_count = points.size;
	MBuf data;
	data.vtx_attr = VtxAttr::PN;
	data.reserve_vertices(points.size + 2); /* +2 for dummy box corners */

	struct Transform transf;
	if (send_points_to_unit_cube(points, data.positions, transf)) {
		printf("Altitude span too large !\n");
		return (-1);
	}

	/* Build kd-tree for filtered point cloud. */
	KdCoords<3> kdcoords{(const float *)data.positions, points.size};
	KdTree<3> kdtree(3, kdcoords, 10);

	/* Build normals using geometry and scanlines*/
	printf("Eval. normal directions    :  ");
	TArray<float> qual(points.size + 2); /* +2 for dummy box corners */
	estim_unoriented_nml(data.positions, points.size, data.normals,
			     qual.data, kdtree, nml_cb);
	printf("\n");

	printf("Eval. normal orientations  :\n");
	TArray<EOrient> oriented(points.size, ENone);

	size_t unset;

	unset =
	    orient_nml_with_scan(points.data, points.size, fls.data, fls.size,
				 qual.data, data.normals, oriented.data);
	printf("    Oriented after scanlines pass       : %.1f %%\n",
	       (points.size - unset) * 100.f / points.size);

	unset = orient_nml_with_z(data.normals, oriented.data, points.size,
				  qual.data);
	printf("    Oriented after positive z pass      : %.1f %%\n",
	       (points.size - unset) * 100.f / points.size);

	float progress = 1.f;
	int pass = 1;
	while (unset && progress > 0.0001f) {
		size_t newly_set =
		    propagate_nml_once(data.positions, points.size, kdtree,
				       qual.data, data.normals, oriented.data);
		progress = (float)newly_set / unset;
		unset -= newly_set;
		printf("    Oriented after propagate pass nbr %d : %.1f %%\n",
		       pass++, (points.size - unset) * 100.f / points.size);
	}

	if (cfg.nml_confidence) {
		/*weight normals according to orientation found and quality */
		for (size_t i = 0; i < points.size; ++i) {
			float weight = (oriented[i] != ENone) ? qual[i] : 0;
			data.normals[i] *= weight;
		}
	} else {
		/* Skip points with no orientation found */
		size_t new_num = 0;
		for (size_t i = 0; i < points.size; ++i) {
			if (oriented[i] != ENone) {
				data.positions[new_num] = data.positions[i];
				data.normals[new_num] = data.normals[i];
				new_num++;
			}
		}
		printf("    %.1f %% vertices discarded due to lack of clear "
		       "orientation.\n",
		       unset * 100.f / points.size);
		mesh.vertex_count = new_num;
	}

	/* Add dummy points for bounding box to perfectly match unit cube */
	data.positions[mesh.vertex_count] = Vec3{0.f, 0.f, 0.f};
	data.normals[mesh.vertex_count] = Vec3{0.f, 0.f, 0.f};
	qual[mesh.vertex_count] = 0;
	mesh.vertex_count++;
	data.positions[mesh.vertex_count] = Vec3{1.f, 1.f, 1.f};
	data.normals[mesh.vertex_count] = Vec3{0.f, 0.f, 0.f};
	qual[mesh.vertex_count] = 0;
	mesh.vertex_count++;

#if DEBUG_NML_CONFIDENCE
	data.add_vtx_attr(VtxAttr::COL);
	for (size_t i = 0; i < mesh.vertex_count; ++i) {
		uint8_t &r = data.colors[i].x;
		uint8_t &g = data.colors[i].y;
		uint8_t &b = data.colors[i].z;
		float bfloat = 1 - 2 * std::abs(qual[i] - 0.5f);
		float rfloat = qual[i] < 0.5f ? 1 - bfloat : 0;
		float gfloat = qual[i] > 0.5f ? 1 - bfloat : 0;
		r = int(255 * rfloat + 0.5f);
		g = int(255 * gfloat + 0.5f);
		b = int(255 * bfloat + 0.5f);
	}
#endif

	write_mesh(mesh, data, cfg, "points.ply");
	write_transform(transf, cfg);

	return (0);
}

/******************************************************************************
 *
 * III. Functions related to surface mesh (post)processing.
 *
 ******************************************************************************/

static void recut_mesh(Mesh &mesh, MBuf &data, const struct Transform &transf)
{
	/* Remove triangles in buffered part */
	/* Note: here we have all offsets 0 so forget about them */
	size_t new_index_count = 0;
	for (size_t i = 0; i < mesh.index_count / 3; ++i) {
		uint32_t i0 = data.indices[3 * i + 0];
		uint32_t i1 = data.indices[3 * i + 1];
		uint32_t i2 = data.indices[3 * i + 2];
		Vec3 p0 = data.positions[i0];
		Vec3 p1 = data.positions[i1];
		Vec3 p2 = data.positions[i2];
		Vec3 bary = (p0 + p1 + p2) * (1 / 3.f);
		if (bary.x < transf.shift.x || bary.x > 1.f - transf.shift.x ||
		    bary.y < transf.shift.y || bary.y > 1.f - transf.shift.y) {
			continue;
		}
		data.indices[new_index_count++] = i0;
		data.indices[new_index_count++] = i1;
		data.indices[new_index_count++] = i2;
	}
	mesh.index_count = new_index_count;
}

static void rescale_and_offset_mesh(Mesh &mesh, MBuf &data,
				    const struct Transform &transf,
				    const struct Cfg &cfg)
{
	/* Inverse transform points + shift relative to base */
	/* TODO : should be eventually removed and use scene
	 *        object placement and scaling instead
	 */
	float invscale = 1. / transf.scale;
	for (size_t i = 0; i < mesh.vertex_count; ++i) {
		data.positions[i] = data.positions[i] - transf.shift;
		data.positions[i] *= invscale;
		data.positions[i].x += cfg.x0 - cfg.x0_base;
		data.positions[i].y += cfg.y0 - cfg.y0_base;
	}
}

static int build_surface_mesh(Mesh &mesh, MBuf &data, const struct Cfg &cfg)
{
	/* Re-use existing output ? */
	char *recon_out =
	    get_filename(cfg.x0, cfg.y0, cfg.out_dir, "poisson.ply");
	FILE *f;
	if ((f = fopen(recon_out, "rb")) != NULL) {
		printf("Using cached data in %s\n", recon_out);
		load_ply(mesh, data, recon_out);
		printf("A total of %d (%.2f M) Tri.\n", mesh.index_count / 3,
		       1e-6 * mesh.index_count / 3);
		free(recon_out);
		fclose(f);
		return 0;
	}

	char *recon_in =
	    get_filename(cfg.x0, cfg.y0, cfg.out_dir, "points.ply");
	const char *format =
	    "poissonrecon --in %s --out %s --scale 1.0 --depth %d "
	    "--pointWeight %d --confidence %d --threads 8 %s";
	unsigned short len =
	    strlen(recon_in) + strlen(recon_out) + strlen(format) + 16;

	char *cmd = (char *)calloc(len, sizeof(*cmd));
	snprintf(cmd, len, format, recon_in, recon_out, cfg.depth, cfg.weight,
		 cfg.nml_confidence ? 1 : 0, cfg.verbose ? "--verbose" : "");
	int ret = system(cmd);

	if (cfg.clean) {
		snprintf(cmd, len, "rm %s", recon_in);
		system(cmd);
	}

	if (!ret) {
		load_ply(mesh, data, recon_out);
		printf(
		    "A total of %d (%.2f M) Tri after poisson reconstruct.\n",
		    mesh.index_count / 3, 1e-6 * mesh.index_count / 3);
		struct Transform transf;
		read_transform(transf, cfg);
		recut_mesh(mesh, data, transf);
		printf(
		    "A total of %d (%.2f M) Tri after buffered boundary cut.\n",
		    mesh.index_count / 3, 1e-6 * mesh.index_count / 3);
		rescale_and_offset_mesh(mesh, data, transf, cfg);
		// clusterize_mesh(mesh, data, 0.001 * 0.1);
		//  optimize_mesh(mesh, data);
		compact_mesh(mesh, data);
		if (cfg.clean) {
			snprintf(cmd, len, "rm %s", recon_out);
			system(cmd);
		} else {
			write_ply(recon_out, mesh, data);
		}
	}

	free(cmd);
	free(recon_out);
	free(recon_in);

	return (ret);
}

static int simplify_mesh(Mesh &mesh, MBuf &data, const struct Cfg &cfg)
{
	float tol = cfg.simp_error;
	mesh.index_count = meshopt_simplify(
	    data.indices, data.indices, mesh.index_count,
	    (const float *)data.positions, mesh.vertex_count, sizeof(Vec3),
	    mesh.index_count * 3 / 4, tol, 0, NULL);

	mesh.vertex_count = meshopt_optimizeVertexFetch(
	    data.positions, data.indices, mesh.index_count, data.positions,
	    mesh.vertex_count, sizeof(Vec3));

	write_mesh(mesh, data, cfg, "meshopt.ply");

	printf("After meshopt simplification : %d (%.2f M) Tri\n",
	       mesh.index_count / 3, 1e-6 * mesh.index_count / 3);

	return (0);
}

#if 0 // Old version
static int improve_mesh_quality(Mesh &mesh, MBuf &data, const struct Cfg &cfg)
{
	char *fin = get_filename(cfg.x0, cfg.y0, cfg.out_dir, "mesh");
	char *fout = get_filename(cfg.x0, cfg.y0, cfg.out_dir, "out.mesh");
	char *fsol = get_filename(cfg.x0, cfg.y0, cfg.out_dir, "out.sol");

	/* TODO Avoid going to disk */
	write_inria(fin, mesh, data);

	const char *format =
	    "mmgs -in %s -out %s -hausd %.8f -hgrad %.8f -nr -v %d %s";
	int len = strlen(fin) + strlen(fout) + strlen(format) + 24;

	char *cmd = (char *)calloc(len, sizeof(*cmd));
	snprintf(cmd, len, format, fin, fout, cfg.hausd * 0.001, cfg.hgrad,
		 cfg.verbose ? 1 : 0, cfg.ani ? "-A" : "");
	int ret = system(cmd);

	snprintf(cmd, len, "rm %s", fin);
	// system(cmd);

	if (!ret) {
		read_inria(mesh, data, fout);
		if (!cfg.clean)
			write_mesh(mesh, data, cfg, "mmgs.ply");
	}

	snprintf(cmd, len, "rm -f %s", fout);
	// system(cmd);
	snprintf(cmd, len, "rm -f %s", fsol);
	// system(cmd);

	free(cmd);
	free(fsol);
	free(fout);
	free(fin);

	return (ret);
}
#else
static int improve_mesh_quality(Mesh &mesh, MBuf &data, const struct Cfg &cfg)
{
	return mesh_remesh(mesh, data, cfg.hausd * 0.001, 5, false, 1);
}
#endif


int quantize_encode_mesh(Mesh &mesh, MBuf &data, const Cfg &cfg)
{
	struct Transform transf;
	if (read_transform(transf, cfg)) {
		printf("Could not read transform\n");
		return (-1);
	}

	uint32_t index_count = mesh.index_count;
	uint32_t vertex_count = mesh.vertex_count;
	TArray<TVec3<uint16_t>> qpos(vertex_count);
	for (size_t i = 0; i < vertex_count; ++i) {
		Vec3 cubepos = data.positions[i];
		cubepos.x -= cfg.x0 - cfg.x0_base;
		cubepos.y -= cfg.y0 - cfg.y0_base;
		cubepos *= transf.scale;
		cubepos += transf.shift;
		qpos[i].x = cubepos.x * ((1 << 16) - 1);
		qpos[i].y = cubepos.y * ((1 << 16) - 1);
		qpos[i].z = cubepos.z * ((1 << 16) - 1);
	}
	uint32_t *indices = data.indices + mesh.index_offset;
	// void *vertices = (float *)(data.positions + mesh.vertex_offset);
	// size_t vertex_size = sizeof(Vec3);
	void *vertices = qpos.data;
	size_t vertex_size = sizeof(TVec3<uint16_t>);
	TArray<uint8_t> vbuf(
	    meshopt_encodeVertexBufferBound(vertex_count, vertex_size));
	vbuf.resize(meshopt_encodeVertexBuffer(&vbuf[0], vbuf.size, vertices,
					       vertex_count, vertex_size));
	TArray<uint8_t> ibuf(
	    meshopt_encodeIndexBufferBound(index_count, vertex_count));
	ibuf.resize(meshopt_encodeIndexBuffer(&ibuf[0], ibuf.size, indices,
					      index_count));

	printf("Sizes : %.1f%% (vertex) %.1fzu%% (index)\n",
	       ((float)vbuf.size / (ibuf.size + vbuf.size)),
	       ((float)ibuf.size / (ibuf.size + vbuf.size)));

	char *fname = get_filename(cfg.x0, cfg.y0, cfg.out_dir, "bin");
	FILE *f = fopen(fname, "wb");
	int ret = (f == NULL) || fwrite(vbuf.data, vbuf.size, 1, f) != 1 ||
			  fwrite(ibuf.data, ibuf.size, 1, f) != 1
		      ? -1
		      : 0;
	fclose(f);
	free(fname);
	return (ret);
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

	printf("\n");
	printf("I. Building oriented point set :\n");
	printf("--------------------------------\n");
	if (build_oriented_point_set(cfg)) {
		printf("No data found.\n");
		return (-1);
	}

	printf("\n");
	printf("II. Building surface mesh from point set :\n");
	printf("------------------------------------------\n");

	if (build_surface_mesh(mesh, data, cfg)) {
		printf("Error in Poisson reconstruction\n");
		return (-1);
	}

	printf("\n");
	printf("III. Postprocessing surface mesh.\n");
	printf("---------------------------------\n");
	if (cfg.simp_error > 0) {
		simplify_mesh(mesh, data, cfg);
	}
	if (cfg.hausd > 0) {
		if (improve_mesh_quality(mesh, data, cfg)) {
			printf("Error in MMGS\n");
			return (-1);
		}
	}
	if (cfg.optimize)
		optimize_mesh(mesh, data);

	/* Save final mesh */
	if (cfg.encode && quantize_encode_mesh(mesh, data, cfg)) {
		printf("Error in mesh encoding\n");
		return (-1);
	}

	write_mesh(mesh, data, cfg, "final.ply");

	printf("\n------ Finished with %04d %04d ------\n", cfg.x0, cfg.y0);
	return (0);
}

/******************************************************************************/
