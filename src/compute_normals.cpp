#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "array.h"
#include "hash_table.h"
#include "kd_tree.h"
#include "math_utils.h"
#include "sys_utils.h"

#include "las_normal.h"
#include "las_read.h"
#include "las_source.h"
#include "mesh.h"
#include "mesh_ply.h"

int inline get_source_idx(const TArray<int> &sources, int source_id)
{
	int n = sources.size;
	while (n--) {
		if (sources[n] == source_id)
			return (n);
	}
	return (-1);
}

inline bool filter_point(const struct LasPoint &p)
{
	return (1); // p.classification == 2);
}

int get_ids_and_counts(TArray<int> &source_ids, TArray<int> &source_counts,
		       const char *fname)
{
	struct LasFileInfo info;
	if (las_read_info(fname, info)) {
		return (-1);
	}
	char *raw_data = las_load_raw_data(fname, info, NULL);
	assert(raw_data);

	struct LasPoint p;
	for (size_t i = 0; i < info.point_num; ++i) {
		const char *raw = raw_data + i * info.point_size;
		p = las_read_point(raw, info.point_format);
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
	return (0);
}

int accumulate_counts(TArray<int> &offsets, const TArray<int> &counts)
{
	assert(offsets.size == counts.size);
	offsets[0] = 0;
	for (size_t i = 1; i < offsets.size; ++i) {
		offsets[i] = offsets[i - 1] + counts[i - 1];
	}
	return offsets[offsets.size - 1] + counts[counts.size - 1];
}

void fill_points(TArray<struct LasPoint> &points, const TArray<int> &ids,
		 const TArray<int> &counts, TArray<int> &offsets,
		 const char *fname)
{
	struct LasFileInfo info;
	las_read_info(fname, info);
	char *raw_data = las_load_raw_data(fname, info, NULL);
	assert(raw_data);
	struct LasPoint p;
	for (size_t i = 0; i < info.point_num; ++i) {
		const char *raw = raw_data + i * info.point_size;
		p = las_read_point(raw, info.point_format);
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
	/* Reset offsets to their correct values */
	for (size_t i = 0; i < offsets.size; ++i) {
		offsets[i] -= counts[i];
	}
}

int read_and_filter_data(size_t &source_num, TArray<struct LasPoint> &points,
			 const char *fname)
{
	TArray<int> ids;
	TArray<int> counts;
	if (get_ids_and_counts(ids, counts, fname))
		return (-1);
	source_num = ids.size;
	if (!source_num)
		return (-1);
	TArray<int> offsets(source_num);
	int point_num = accumulate_counts(offsets, counts);
	points.clear();
	points.resize(point_num);
	fill_points(points, ids, counts, offsets, fname);

	return (0);
}

int set_positions(const TArray<struct LasPoint> &points,
		  const LasFileInfo &info, Vec3 *pos)
{
	for (size_t i = 0; i < points.size; ++i) {
		pos[i].x = points[i].x * info.scale[0];
		pos[i].y = points[i].y * info.scale[1];
		pos[i].z = points[i].z * info.scale[2];
	}
	return (0);
}

static void nml_cb(float progress)
{
	printf("\rBuilding unoriented normals : %.0f %%", progress);
	fflush(stdout);
}

int main(int argc, char **argv)
{
	/* Process command line arguments */
	if (argc < 2) {
		printf("Syntax : compute_normals filename\n");
		return (-1);
	}
	const char *fname = argv[1];

	/* Read info */
	struct LasFileInfo info;
	if (las_read_info(fname, info)) {
		printf("Error reading %s\n", fname);
		return (-1);
	}

	/* Read data */
	size_t source_num;
	TArray<struct LasPoint> points;
	if (read_and_filter_data(source_num, points, fname)) {
		printf("Data contains no point.\n");
		return (-1);
	}
	size_t point_num = points.size;

	printf("Total sources : %zu\n", source_num);
	printf("Total filtered points : %zu\n", point_num);

	/* Derive fligh lines azimut from data. Used later for normals
	 */
	printf("Reconstructing flightlines.\n");
	TArray<struct SourceStat> stats(source_num);
	las_stat_sources(points, stats);
	TArray<struct SourceFlightLine> fls(source_num);
	las_approx_flight_lines(points, info.scale, stats, fls);

	printf("Rescaling positions to unit cube.\n");
	Mesh mesh;
	MBuf data;
	mesh.vertex_count = point_num;
	data.vtx_attr = VtxAttr::PN;
	data.reserve_vertices(point_num);
	set_positions(points, info, data.positions);

	/* Build kd-tree for filtered point cloud.
	 */
	printf("Building kd-tree.\n");
	const float *pos = (const float *)data.positions;
	KdTree<3> tree(3, KdCoords<3>{pos, point_num}, 10);

	/* Build normals using geometry and scanlines*/
	TArray<float> qual(point_num);
	printf("Building unoriented normals :");
	estim_unoriented_nml(data.positions, point_num, data.normals, qual.data,
			     tree, nml_cb);
	printf("\n");
	printf("Orienting normals :\n");
	TArray<EOrient> oriented(point_num, ENone);
	size_t unsettled;
	if (argc < 3) {
		unsettled = orient_nml_with_scan(
		    points.data, point_num, fls.data, source_num, qual.data,
		    data.normals, oriented.data);
		printf(
		    "    Oriented after scanlines orientation pass : %.0f %%\n",
		    (point_num - unsettled) * 100.0 / point_num);
	}
	unsettled = orient_nml_with_z(data.normals, oriented.data, point_num,
				      qual.data);
	printf("    Oriented after positive z pass : %.0f %%\n",
	       (point_num - unsettled) * 100.0 / point_num);
	float progress = 1;
	int pass = 1;
	while (unsettled && (progress > 0.0001 || pass < 10) && pass < 50) {
		size_t newly_settled = propagate_nml_once(
		    data.positions, point_num, tree, qual.data, data.normals,
		    oriented.data, 20, 0.7);
		progress = (float)newly_settled / unsettled;
		unsettled -= newly_settled;
		printf("\r    Oriented after propagate pass %d : %.1f %%\n",
		       pass++, (point_num - unsettled) * 100.0 / point_num);
		fflush(stdout);
	}
	printf("\n");

#if 1
	/*weight normals according to orientation found and quality */
	for (size_t i = 0; i < points.size; ++i) {
		float weight = (oriented[i] != ENone) ? qual[i] : 0;
		data.normals[i] *= weight;
	}
#else
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
#endif
#if 1
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

	write_ply("normal_test.ply", mesh, data);

	return (0);
}
