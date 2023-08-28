#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "meshoptimizer/src/meshoptimizer.h"
#include "mmg/mmgs/libmmgs.h"

#include "array.h"
#include "chrono.h"
#include "hash_table.h"
#include "math_utils.h"
#include "sys_utils.h"

#include "mesh.h"
#include "mesh_adjacency.h"
#include "mesh_inria.h"
#include "mesh_ply.h"
#include "mesh_remesh.h"
#include "mesh_utils.h"
#include "vertex_table.h"

#include "copc.h"
#include "las_normal.h"
#include "las_read.h"
#include "las_source.h"

//#include "swiss_raster.h"

/* Size of tile boundary buffer in cm */
#define BDY_BUFFER 10000
#define BDY_BUFFER_ADD 48000

/* Save normal quality as PLY point color attribute */
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
	float weight;
	float hausd;
	float hgrad;
	int clean;
	bool verbose;
	bool optimize;
	bool encode;
	bool nml_confidence;
};

struct Timings {
	unsigned int read_and_filter = 0;
	unsigned int estim_nml = 0;
	unsigned int poisson_recon = 0;
	unsigned int mmgs = 0;
	unsigned int total = 0;
};

static int process_args(int argc, const char **argv, struct Cfg &cfg)
{
	if (argc < 3) {
		printf("Syntax : french_lidar x y [...]\n");
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
	cfg.weight = (argc >= 7) ? atof(argv[6]) : 4;
	cfg.hausd = (argc >= 8) ? atof(argv[7]) : 0.1f;
	cfg.hgrad = (argc >= 9) ? atof(argv[8]) : 5.f;
	cfg.clean = (argc >= 10) ? atoi(argv[9]) : 1;
	cfg.verbose = (argc >= 11) ? atoi(argv[10]) : 1;
	cfg.optimize = (argc >= 12) ? atoi(argv[11]) : 1;
	cfg.encode = (argc >= 13) ? atoi(argv[12]) : 0;
	cfg.nml_confidence = (argc >= 14) ? atoi(argv[13]) : 1;

	return (0);
}

static void print_cfg(const struct Cfg &cfg)
{
	printf("\n");
	printf("Configuration :\n");
	printf("---------------\n");
	printf("Tile coords : %d %d\n", cfg.x0, cfg.y0);
	printf("Data  dir   : %s\n", cfg.base_dir);
	printf("Output dir  : %s\n", cfg.out_dir);
	printf("Verbosity   : %d\n", cfg.verbose ? 1 : 0);
	printf("Simplify    : hausd=%f , hgrad=%f\n", cfg.hausd, cfg.hgrad);
}

static void print_timings(const Timings &tt)
{
	unsigned div = 1000000;
	unsigned int other = tt.total - tt.read_and_filter - tt.estim_nml -
			     tt.poisson_recon - tt.mmgs;
	printf("\n");
	printf("Timings :\n");
	printf("---------\n");
	printf("Read data   : ");
	printf("%3d s (%.1f%%)\n", tt.read_and_filter / div,
	       100.f * tt.read_and_filter / tt.total);
	printf("Estim nml   : ");
	printf("%3d s (%.1f%%)\n", tt.estim_nml / div,
	       100.f * tt.estim_nml / tt.total);
	printf("Poisson     : ");
	printf("%3d s (%.1f%%)\n", tt.poisson_recon / div,
	       100.f * tt.poisson_recon / tt.total);
	printf("MMGS        : ");
	printf("%3d s (%.1f%%)\n", tt.mmgs / div, 100.f * tt.mmgs / tt.total);
	printf("Other       : ");
	printf("%3d s (%.1f%%)\n", other / div, 100.f * other / tt.total);
	printf("Total       : ");
	printf("%3d s \n", tt.total / div);
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
	snprintf(fname + base_len, 32, "%04d_%04d.%s", x, y, ext);

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

static bool filter_las_point(const LasPoint &p, const LasFileInfo &info,
			     const TAabb<double> box)
{
	double pos[3];
	pos[0] = p.x * info.scale[0] + info.offset[0];
	pos[1] = p.y * info.scale[1] + info.offset[1];
	pos[2] = p.z * info.scale[2] + info.offset[2];

	/* Bbox filter */
	if ((pos[0] < box.min[0]) || (pos[0] > box.max[0]) ||
	    (pos[1] < box.min[1]) || (pos[1] > box.max[1])) {
		return false;
	}

	/* Altitude filter */
	if (pos[2] < 2200) {
		return (p.classification == 2 || p.classification == 9);
	} else {
		return p.classification <= 9;
	}
}

#if 0
static bool filter_raster_point(const struct LasPoint &p,
				const TArray<uint32_t> &density, int N)
{
	/* 100m boundary buffer size */
	const int bd = BDY_BUFFER;
	if ((p.x < -bd) || (p.x > 100000 + bd) || (p.y < -bd) ||
	    (p.y > 100000 + bd)) {
		return false;
	}
	/* Density check */
	const float width = 100000 + 2 * bd;
	const float invh = N / width;
	int pix_x = ((float)p.x + bd) * invh - 0.5f;
	pix_x = pix_x < 0 ? 0 : pix_x;
	pix_x = pix_x >= N ? N - 1 : pix_x;
	int pix_y = ((float)p.y + bd) * invh - 0.5f;
	pix_y = pix_y < 0 ? 0 : pix_y;
	pix_y = pix_y >= N ? N - 1 : pix_y;
	int pix = pix_y * N + pix_x;
	return (density[pix] == 0);
}
#endif

#if 0
static bool filter_bdy_raster_point(const struct LasPoint &p, int N)
{
	/* 500m to 100m boundary buffer size */
	const int bd_in = BDY_BUFFER;
	const int bd_out = BDY_BUFFER_ADD;
	if ((p.x <= -bd_out) || (p.x >= 100000 + bd_out) || (p.y <= -bd_out) ||
	    (p.y >= 100000 + bd_out)) {
		return false;
	}
	if ((p.x >= -bd_in) && (p.x <= 100000 + bd_in) && (p.y >= -bd_in) &&
	    (p.y <= 100000 + bd_in)) {
		return false;
	}
	return (true);
}
#endif

static TAabb<double> las_bbox(int x0, int y0)
{
	TAabb<double> box;
	box.min.x = 1000 * x0 - 100;
	box.min.y = 1000 * y0 - 1100;
	box.min.z = -1000;
	box.max.x = 1000 * x0 + 1100;
	box.max.y = 1000 * y0 + 100;
	box.max.z = 9000;
	return box;
}

static uint32_t filter_and_add_points(TArray<LasPoint> &points, const char *src,
				      uint32_t src_count,
				      const LasFileInfo info,
				      const TAabb<double> box)
{
	size_t init_point_count = points.size;
	size_t point_count = init_point_count;
	points.resize(point_count + src_count); // Over estimate, shrink later
	for (uint32_t i = 0; i < src_count; ++i) {
		LasPoint &p = points[point_count];
		p = las_read_point(src, info.point_format);
		if (filter_las_point(p, info, box))
			point_count++;
		src += info.point_size;
	}
	points.resize(point_count);

	return (point_count - init_point_count);
}

static size_t read_and_filter_las_data(TArray<struct LasPoint> &points,
				       double offset[3], double scale[3],
				       const struct Cfg &cfg)
{
	TAabb<double> box = las_bbox(cfg.x0, cfg.y0);
	TArray<char> buf;
	for (int i = 0; i < 9; ++i) {
		if (cfg.verbose) {
			printf("\rReading point cloud data [%2d]", i);
			fflush(stdout);
		}
		int dx = (i % 3) - 1;
		int dy = (i / 3) - 1;
		int x = cfg.x0 + dx;
		int y = cfg.y0 + dy;
		char *fname = get_filename(x, y, cfg.base_dir, "copc.laz");
		LasFileInfo info;
		if (las_read_info(fname, info)) {
			free(fname);
			continue;
		}
		assert(info.copc);
		struct CopcReader *copc = copc_init(fname);
		if (!copc) {
			free(fname);
			continue;
		}
		uint32_t cell_count = copc_set_target_bbox(copc, box);
		// printf("Cell count : %d\n", cell_count);
		for (uint32_t k = 0; k < cell_count; ++k) {
			int cell_points = copc_cell_point_count(copc, k);
			assert(cell_points);
			buf.reserve(cell_points * info.point_size);
			copc_read_cell(copc, k, &buf[0]);
			filter_and_add_points(points, &buf[0], cell_points,
					      info, box);
			// printf("Cell points : %d, New points : %d\n",
			//        cell_points, new_points);
		}
		copc_fini(copc);
		free(fname);
	}
	if (cfg.verbose)
		printf("\n");
	return points.size;
}

static size_t analyze_source_flight_lines(TArray<struct LasPoint> &points,
					  TArray<struct SourceFlightLine> &fls,
					  const struct Cfg &cfg)
{
	TArray<int> sources;
	for (size_t i = 0; i < points.size; ++i) {
		LasPoint &p = points[i];
		int idx = get_source_idx(sources, p.source_id);
		if (idx >= 0) {
			p.source_idx = idx;
		} else {
			p.source_idx = sources.size;
			sources.push_back(p.source_id);
		}
	}
	/* Derive fligh lines azimut from data. Used later for normals
	 */
	if (cfg.verbose)
		printf("Reconstructing flightlines :");
	TArray<struct SourceStat> stats(sources.size);
	las_stat_sources(points, stats);
	fls.resize(sources.size);
	double scale[3] = {1., 1., 1.};
	int valid = las_approx_flight_lines(points, scale, stats, fls);
	if (cfg.verbose)
		printf(" found %d valid out of %zu sources\n", valid,
		       sources.size);
	return sources.size;
}

#if 0
static size_t get_raster_point_count(int x0, int y0, const char *base_dir,
				     const TArray<uint32_t> &density, int N)
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
				/* No data check */
				if (altitudes[nx * py + px] < -500)
					continue;
				p.z = (int)(100 * altitudes[nx * py + px]);
				rebase_las_point(p, dx, dy);
				if (filter_raster_point(p, density, N))
					++raster_point_count;
			}
		}
	}
	return (raster_point_count);
}

static size_t fill_raster_points(TArray<struct LasPoint> &points, size_t offset,
				 int dummy_source_idx, int x0, int y0,
				 const char *base_dir,
				 const TArray<uint32_t> &density, int N)
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
				/* No data check */
				if (altitudes[nx * py + px] < -500)
					continue;
				p.z = (int)(100 * altitudes[nx * py + px]);
				p.source_idx = dummy_source_idx;
				rebase_las_point(p, dx, dy);
				if (filter_raster_point(p, density, N))
					points[point_idx++] = p;
			}
		}
	}
	return (point_idx - offset);
}

static size_t fill_bdy_raster_points(TArray<struct LasPoint> &points,
				     size_t offset, int dummy_source_idx,
				     int x0, int y0, const char *base_dir,
				     int min_z, int max_z, int N)
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
				/* No data check */
				if (altitudes[nx * py + px] < -500)
					continue;
				p.z = (int)(100 * altitudes[nx * py + px]);
				p.z = p.z > max_z ? max_z : p.z;
				p.z = p.z < min_z ? min_z : p.z;
				p.source_idx = dummy_source_idx;
				rebase_las_point(p, dx, dy);
				if (filter_bdy_raster_point(p, N))
					points[point_idx++] = p;
			}
		}
	}
	return (point_idx - offset);
}
#endif

#if 0
static void
build_las_point_density_matrix(const TArray<struct LasPoint> &points,
			       TArray<uint32_t> &density, int N)
{
	density.clear();
	density.resize(N * N);
	memset(&density[0], 0, N * N * sizeof(uint32_t));
	float width = 100000 + 2 * BDY_BUFFER;
	float invh = N / width;
	for (size_t i = 0; i < points.size; ++i) {
		int pix_x = ((float)points[i].x + BDY_BUFFER) * invh - 0.5f;
		pix_x = pix_x < 0 ? 0 : pix_x;
		pix_x = pix_x >= N ? N - 1 : pix_x;
		int pix_y = ((float)points[i].y + BDY_BUFFER) * invh - 0.5f;
		pix_y = pix_y < 0 ? 0 : pix_y;
		pix_y = pix_y >= N ? N - 1 : pix_y;
		int pix = pix_y * N + pix_x;
		density[pix] += 1;
	}
}
#endif

#if 0
static size_t read_and_filter_raster_data(TArray<struct LasPoint> &points,
					  TArray<struct SourceFlightLine> &fls,
					  const struct Cfg &cfg)
{
	size_t raster_point_count = 0;

	const uint32_t N = (100000 + 2 * BDY_BUFFER) / 200;
	TArray<uint32_t> density;
	build_las_point_density_matrix(points, density, N);

	/* Get hole filling raster points */
	uint32_t empty_pix = 0;
	for (size_t pix = 0; pix < N * N; pix++) {
		if (!density[pix])
			empty_pix++;
	}
	if (empty_pix) {
		raster_point_count = get_raster_point_count(
		    cfg.x0, cfg.y0, cfg.base_dir, density, N);
	}
	if (raster_point_count) {
		size_t offset = points.size;
		points.reserve(points.size + raster_point_count);
		points.resize(points.size + raster_point_count);
		size_t dummy_source_idx = fls.size;
		fls.resize(fls.size + 1);
		fls[fls.size - 1].is_valid = false;
		fill_raster_points(points, offset, dummy_source_idx, cfg.x0,
				   cfg.y0, cfg.base_dir, density, N);
	}

	/* Enlarge boundary buffer ? */
	int min_z = INT_MAX, max_z = -INT_MAX;
	for (size_t i = 0; i < points.size; ++i) {
		min_z = MIN(min_z, points[i].z);
		max_z = MAX(max_z, points[i].z);
	}
	float span = (max_z - min_z) * 0.01f;
	if (span >= 1250.f) {
		size_t offset = points.size;
		size_t n_out = (100000 + 2 * BDY_BUFFER_ADD) / 200;
		size_t n_in = (100000 + 2 * BDY_BUFFER) / 200;
		size_t max_bdy_raster_point =
		    (n_out + n_in + 1) * (n_out - n_in + 1);
		points.reserve(points.size + max_bdy_raster_point);
		points.resize(points.size + max_bdy_raster_point);
		if (!raster_point_count) {
			fls.resize(fls.size + 1);
			fls[fls.size - 1].is_valid = false;
		}
		size_t bdy_raster_point_count = fill_bdy_raster_points(
		    points, offset, fls.size - 1, cfg.x0, cfg.y0, cfg.base_dir,
		    min_z, max_z, N);
		assert(bdy_raster_point_count <= max_bdy_raster_point);
		points.resize(offset + bdy_raster_point_count);
		raster_point_count += bdy_raster_point_count;
	}

	return raster_point_count;
}
#endif

static int send_points_to_unit_cube(const TArray<struct LasPoint> &points,
				    Vec3 *pos, struct Transform &t,
				    const Cfg cfg)
{
	int min_z = INT_MAX, max_z = -INT_MAX;
	for (size_t i = 0; i < points.size; ++i) {
		min_z = MIN(min_z, points[i].z);
		max_z = MAX(max_z, points[i].z);
	}
	float span = (max_z - min_z) * 0.01f;
	printf("(Altitude span : %.0f meters)\n", span);
	float n;
	if (span < 1250.f) {
		n = 6;
	} else if (span < 1875.f) {
		n = 2;
	} else {
		return (-1);
	}
	t.scale = (float)n / (n + 2);
	t.shift.x = t.shift.y = 1.f / (n + 2);
	float mean = 0.5 * (min_z + max_z) * 1e-5 * t.scale;
	/* Round shift.z so that octree boxes match vertically to
	 * neighboors */
	t.shift.z = 0.5 - round(16 * mean) * 0.0625;
	assert(min_z * 1e-5 * t.scale + t.shift.z >= 0);
	assert(max_z * 1e-5 * t.scale + t.shift.z <= 1);
	for (size_t i = 0; i < points.size; ++i) {
		double scal = 1e-5 * t.scale;
		pos[i].x = (points[i].x - 100000 * cfg.x0) * scal + t.shift.x;
		pos[i].y =
		    (points[i].y - 100000 * (cfg.y0 - 1)) * scal + t.shift.y;
		pos[i].z = points[i].z * scal + t.shift.z;
		// printf("%lf %lf %lf\n", pos[i].x, pos[i].y, pos[i].z);
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
static int build_oriented_point_set(const struct Cfg &cfg, Timings &tt)
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
			return (0);
		}
	}

	/* Read data */
	Timer chrono;
	chrono.start();
	TArray<struct LasPoint> points;
	TArray<struct SourceFlightLine> fls;
	double offset[3];
	double scale[3];
	size_t las_num = read_and_filter_las_data(points, offset, scale, cfg);
	tt.read_and_filter = chrono.stop();

	/* Analyze flight lines */
	analyze_source_flight_lines(points, fls, cfg);

	/* Read (and filter) raster data */
	// size_t ras_num = read_and_filter_raster_data(points, fls, cfg);
	size_t ras_num = 0;

	printf("Total Lidar points used    : %zu\n", las_num);
	// printf("Total Raster points used   : %zu\n", ras_num);

	if (!las_num && !ras_num) {
		printf("No data available for this tile.\n");
		return (-1);
	}

	printf("Set positions & transform  : ");
	/* Rescale and offset positions into buffer */
	Mesh mesh;
	mesh.vertex_count = points.size;
	MBuf data;
	data.vtx_attr = VtxAttr::PN;
	data.reserve_vertices(points.size + 2); /* +2 for dummy box corners */

	struct Transform transf;
	if (send_points_to_unit_cube(points, data.positions, transf, cfg)) {
		printf("Altitude span too large !\n");
		return (-1);
	}

	/* Build kd-tree for filtered point cloud. */
	KdCoords<3> kdcoords{(const float *)data.positions, points.size};
	KdTree<3> kdtree(3, kdcoords, 10);

	/* Build normals using geometry and scanlines*/
	chrono.start();
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
	printf("    Oriented after scanlines pass        : %.1f %%\n",
	       (points.size - unset) * 100.f / points.size);

	unset = orient_nml_with_z(data.normals, oriented.data, points.size,
				  qual.data);
	printf("    Oriented after positive z pass       : %.1f %%\n",
	       (points.size - unset) * 100.f / points.size);

	int pass = 1;
	float progress = 1.f;
	while (unset && progress && pass < 100) {
		size_t newly_set;
		newly_set =
		    propagate_nml_once(data.positions, points.size, kdtree,
				       qual.data, data.normals, oriented.data);
		progress = (float)newly_set / unset;
		unset -= newly_set;
		printf("\r    Oriented after propagate pass nbr %2d : "
		       "%.1f %%",
		       pass++, (points.size - unset) * 100.f / points.size);
		fflush(stdout);
	}
	printf("\n");

	if (cfg.nml_confidence) {
		/*weight normals according to orientation found and
		 * quality */
		for (size_t i = 0; i < points.size; ++i) {
			float weight = (oriented[i] >= EOriented) ? qual[i] : 0;
			data.normals[i] *= weight;
		}
		printf("There are %zu points remaining without clear "
		       "normal "
		       "orientation.\nTheir normal was set to zero.\n",
		       unset);
	} else {
		/* Skip points with no orientation found */
		size_t new_num = 0;
		for (size_t i = 0; i < points.size; ++i) {
			if (oriented[i] >= EOriented) {
				data.positions[new_num] = data.positions[i];
				data.normals[new_num] = data.normals[i];
				new_num++;
			}
		}
		printf("There were %zu vertices discarded due to lack "
		       "of clear "
		       "orientation.\n",
		       unset);
		mesh.vertex_count = new_num;
	}
	tt.estim_nml = chrono.stop();

	/* Add dummy points for bounding box to perfectly match unit
	 * cube */
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
	}
}

static int build_surface_mesh(const struct Cfg &cfg, Timings &tt)
{
	/* Re-use existing output ? */
	char *recon_out =
	    get_filename(cfg.x0, cfg.y0, cfg.out_dir, "poisson.ply");
	FILE *f;
	if ((f = fopen(recon_out, "rb")) != NULL) {
		printf("Using cached data in %s\n", recon_out);
		free(recon_out);
		fclose(f);
		return 0;
	}

	Timer chrono;
	chrono.start();
	char *recon_in =
	    get_filename(cfg.x0, cfg.y0, cfg.out_dir, "points.ply");
	const char *format =
	    "poissonrecon --in %s --out %s --scale 1.0 --depth %d "
	    "--pointWeight %.1f --confidence %d --threads 8 %s";
	unsigned short len =
	    strlen(recon_in) + strlen(recon_out) + strlen(format) + 32;

	char *cmd = (char *)calloc(len, sizeof(*cmd));
	snprintf(cmd, len, format, recon_in, recon_out, cfg.depth, cfg.weight,
		 cfg.nml_confidence ? 1 : 0, cfg.verbose ? "--verbose" : "");
	int ret = system(cmd);

	if (cfg.clean >= 2) {
		snprintf(cmd, len, "rm -f %s", recon_in);
		system(cmd);
	}

	free(cmd);
	free(recon_out);
	free(recon_in);

	tt.poisson_recon = chrono.stop();
	return (ret);
}

#if 0
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
#endif

uint32_t select_principal_connected_component(Mesh &mesh, MBuf &data)
{
	struct EAdj eadj;
	fill_edge_adjacency(mesh, data, eadj);

	TArray<uint32_t> triadj;
	fill_tri_adjacency(mesh, data, eadj, triadj);

	TArray<uint32_t> cc;
	uint32_t num_cc = find_connected_components(triadj, cc);

	if (num_cc == 1)
		return (num_cc);

	size_t tri_count = mesh.index_count / 3;
	assert(cc.size == tri_count);

	TArray<uint32_t> counts(num_cc, 0);
	for (size_t i = 0; i < tri_count; ++i) {
		counts[cc[i]]++;
	}
	uint32_t cc_max = 0;
	for (size_t i = 1; i < num_cc; ++i) {
		cc_max = counts[i] > counts[cc_max] ? i : cc_max;
	}
	size_t new_index_count = 0;
	uint32_t *indices = data.indices + mesh.index_offset;
	for (size_t i = 0; i < tri_count; ++i) {
		if (cc[i] == cc_max) {
			indices[new_index_count++] = indices[3 * i + 0];
			indices[new_index_count++] = indices[3 * i + 1];
			indices[new_index_count++] = indices[3 * i + 2];
		}
	}
	mesh.index_count = new_index_count;

	return (num_cc);
}
size_t fix_boundary_vertices(const Mesh &mesh, MBuf &data)
{
	TArray<bool> is_bd(mesh.vertex_count, false);
	TArray<uint32_t> counts(mesh.vertex_count, 0);
	TArray<uint32_t> offsets(mesh.vertex_count);
	TArray<uint32_t> edges(mesh.index_count);

	/* Fill counts */
	const uint32_t *idx = data.indices + mesh.index_offset;
	for (size_t i = 0; i < mesh.index_count; ++i) {
		counts[idx[i]] += 1;
	}
	/* Fill offsets */
	uint32_t offset = 0;
	for (size_t i = 0; i < mesh.vertex_count; ++i) {
		offsets[i] = offset;
		offset += counts[i];
	}
	/* Fill edges */
	for (size_t i = 0; i < mesh.index_count / 3; ++i) {
		uint32_t i0 = idx[3 * i + 0];
		uint32_t i1 = idx[3 * i + 1];
		uint32_t i2 = idx[3 * i + 2];
		edges[offsets[i0]++] = i1;
		edges[offsets[i1]++] = i2;
		edges[offsets[i2]++] = i0;
	}
	for (size_t i = 0; i < mesh.vertex_count; ++i) {
		offsets[i] -= counts[i];
	}
	for (size_t i0 = 0; i0 < mesh.vertex_count; ++i0) {
		if (is_bd[i0])
			continue;
		for (size_t j = 0; j < counts[i0]; ++j) {
			uint32_t i1 = edges[offsets[i0] + j];
			bool bd_edge = true;
			for (size_t k = 0; k < counts[i1]; ++k) {
				uint32_t i0b = edges[offsets[i1] + k];
				if (i0 == i0b) {
					bd_edge = false;
					break;
				}
			}
			if (bd_edge) {
				is_bd[i0] = is_bd[i1] = true;
			}
		}
	}
	size_t bd_count = 0;
	size_t interior_bd = 0;
	float max_bd_dist = 0;
	Vec3 *pos = data.positions + mesh.vertex_offset;
	for (size_t i = 0; i < mesh.vertex_count; ++i) {
		if (is_bd[i]) {
			bd_count += 1;
			float errx = std::abs(pos[i].x - roundf(pos[i].x));
			float erry = std::abs(pos[i].y - roundf(pos[i].y));
			if (MIN(errx, erry) > 0.01) {
				interior_bd += 1;
				continue;
			}
			if (errx <= erry) {
				pos[i].x = roundf(pos[i].x);
				max_bd_dist = std::max(max_bd_dist, errx);
			} else {
				pos[i].y = roundf(pos[i].y);
				max_bd_dist = std::max(max_bd_dist, erry);
			}
		}
	}
	printf("Max_bd_dist : %f Interior bd count : %zu\n", max_bd_dist,
	       interior_bd);
	return bd_count;
}

static int improve_mesh_quality(Mesh &mesh, MBuf &data, const struct Cfg &cfg)
{
	/* Probably safer and might help MMGS */
	compact_mesh(mesh, data);

	/* MMGS treatment */
	MMG5_pMesh mm = NULL;
	MMG5_pSol ss = NULL;
	MMGS_Init_mesh(MMG5_ARG_start, MMG5_ARG_ppMesh, &mm, MMG5_ARG_ppMet,
		       &ss, MMG5_ARG_end);

	MMGS_Init_parameters(mm);

	MMGS_Set_dparameter(mm, ss, MMGS_DPARAM_hausd, cfg.hausd * 0.001);
	MMGS_Set_dparameter(mm, ss, MMGS_DPARAM_hgrad, cfg.hgrad);
	MMGS_Set_iparameter(mm, ss, MMGS_IPARAM_angle, 0);
	MMGS_Set_iparameter(mm, ss, MMGS_IPARAM_verbose, cfg.verbose ? 1 : -1);

	load_mesh_to_mmg(mesh, data, mm, ss);

	const Vec3 *pos = data.positions + mesh.vertex_offset;
	for (size_t i = 0; i < mesh.vertex_count; ++i) {
		bool req = false;
		float eps = 0.005;

		req |= (pos[i].x == 0 && pos[i].y <= eps);
		req |= (pos[i].x == 0 && pos[i].y >= 1 - eps);

		req |= (pos[i].x == 1 && pos[i].y <= eps);
		req |= (pos[i].x == 1 && pos[i].y >= 1 - eps);

		req |= (pos[i].y == 0 && pos[i].x <= eps);
		req |= (pos[i].y == 0 && pos[i].x >= 1 - eps);

		req |= (pos[i].y == 1 && pos[i].x <= eps);
		req |= (pos[i].y == 1 && pos[i].x >= 1 - eps);

		if (req)
			MMGS_Set_requiredVertex(mm, i + 1);
	}

	MMGS_mmgslib(mm, ss);
	int np, nt, na;
	MMGS_Get_meshSize(mm, &np, &nt, &na);

	unload_mesh_from_mmg(mesh, data, mm, ss);

	MMGS_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh, &mm, MMG5_ARG_ppMet, &ss,
		      MMG5_ARG_end);

	/* MMGS will move some vertices out of [0.1] (close to corners),
	 * we reproject them to the cube afterwards.
	 */
	fix_boundary_vertices(mesh, data);

	return 0;
}

int write_encoded_mesh(const Mesh &mesh, const MBuf &data, const Cfg &cfg,
		       const char *ext)
{
	uint32_t index_count = mesh.index_count;
	uint32_t vertex_count = mesh.vertex_count;
	TArray<TVec3<uint16_t>> qpos(vertex_count);
	for (size_t i = 0; i < vertex_count; ++i) {
		Vec3 pos = data.positions[i + mesh.vertex_offset];
		qpos[i].x = pos.x * (1 << 15) + (1 << 14);
		qpos[i].y = pos.y * (1 << 15) + (1 << 14);
		qpos[i].z = pos.z * (1 << 14); /* TODO : scale in z */
	}
	uint32_t *indices = data.indices + mesh.index_offset;
	// void *vertices = qpos.data;
	// size_t vertex_size = sizeof(TVec3<uint16_t>);
	void *vertices = data.positions + mesh.vertex_offset;
	size_t vertex_size = sizeof(Vec3);
	TArray<uint8_t> vbuf(
	    meshopt_encodeVertexBufferBound(vertex_count, vertex_size));
	vbuf.resize(meshopt_encodeVertexBuffer(&vbuf[0], vbuf.size, vertices,
					       vertex_count, vertex_size));
	printf("Bytes per vertex : %.1f\n", (float)vbuf.size / vertex_count);
	TArray<uint8_t> ibuf(
	    meshopt_encodeIndexBufferBound(index_count, vertex_count));
	ibuf.resize(meshopt_encodeIndexBuffer(&ibuf[0], ibuf.size, indices,
					      index_count));
	printf("Index bytes per triangle : %.1f\n",
	       3 * (float)ibuf.size / index_count);

	char *fname = get_filename(cfg.x0, cfg.y0, cfg.out_dir, ext);
	FILE *f = fopen(fname, "wb");
	int ret = (f == NULL) || fwrite(vbuf.data, vbuf.size, 1, f) != 1 ||
			  fwrite(ibuf.data, ibuf.size, 1, f) != 1
		      ? -1
		      : 0;
	fclose(f);
	free(fname);
	return (ret);
}

int postprocess_surface_mesh(const Cfg &cfg, Timings &tt)
{
	char *recon_out =
	    get_filename(cfg.x0, cfg.y0, cfg.out_dir, "poisson.ply");

	/* Load from Poisson Recon */
	Mesh mesh;
	MBuf data;
	load_ply(mesh, data, recon_out);
	printf("A total of %d (%.2f M) Tri after poisson reconstruct.\n",
	       mesh.index_count / 3, 1e-6 * mesh.index_count / 3);

	/* Remove spurious connected components if any */
	int num_cc = select_principal_connected_component(mesh, data);
	if (num_cc != 1) {
		printf("Removed %d connected components.\n", num_cc - 1);
		printf("A total of %d (%.2f M) Tri after spurious cc "
		       "removal.\n",
		       mesh.index_count / 3, 1e-6 * mesh.index_count / 3);
	}

	/* Recut mesh to 1km boundary */
	struct Transform transf;
	read_transform(transf, cfg);
	recut_mesh(mesh, data, transf);
	printf("A total of %d (%.2f M) Tri after buffered boundary "
	       "recut.\n",
	       mesh.index_count / 3, 1e-6 * mesh.index_count / 3);

	/* Rescale and offset (scale is now 1 = 1km for x, y and z) */
	rescale_and_offset_mesh(mesh, data, transf, cfg);

	if (cfg.hausd > 0) {
		Timer chrono;
		chrono.start();
		if (improve_mesh_quality(mesh, data, cfg)) {
			printf("Error in MMGS\n");
			return (-1);
		}
		tt.mmgs = chrono.stop();
	}

	if (cfg.optimize) {
		/* Necessary to compact ? */
		compact_mesh(mesh, data);
		optimize_mesh(mesh, data);
	}

	/* Save final mesh */
	if (cfg.encode) {
		write_encoded_mesh(mesh, data, cfg, "final.bin");
	} else {
		write_mesh(mesh, data, cfg, "final.ply");
	}

	if (cfg.clean) {
		unsigned short len = strlen(recon_out) + 8;
		char *cmd = (char *)calloc(len, sizeof(*cmd));
		snprintf(cmd, len, "rm -f %.256s", recon_out);
		system(cmd);
		free(cmd);
	}

	free(recon_out);

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
	struct Timings tt;
	Timer chrono;

	chrono.start();

	/* Process command line arguments */
	if (process_args(argc, (const char **)argv, cfg)) {
		return (-1);
	}

	printf("\n------ Start of french_lidar for %04d %04d ------\n", cfg.x0,
	       cfg.y0);

	if (cfg.verbose) {
		print_cfg(cfg);
	}

	printf("\n");
	printf("I. Building oriented point set :\n");
	printf("--------------------------------\n");
	if (build_oriented_point_set(cfg, tt)) {
		printf("Error in building oriented point set.\n");
		return (-1);
	}

	printf("\n");
	printf("II. Building surface mesh from point set :\n");
	printf("------------------------------------------\n");

	if (build_surface_mesh(cfg, tt)) {
		printf("Error in Poisson reconstruction\n");
		return (-1);
	}

	printf("\n");
	printf("III. Postprocessing surface mesh.\n");
	printf("---------------------------------\n");

	if (postprocess_surface_mesh(cfg, tt)) {
		printf("Error in Poisson reconstruction\n");
		return (-1);
	}

	tt.total = chrono.stop();

	if (cfg.verbose) {
		print_timings(tt);
	}

	printf("\n------ End of french_lidar for %04d %04d ------\n", cfg.x0,
	       cfg.y0);

	return (0);
}

/******************************************************************************/
