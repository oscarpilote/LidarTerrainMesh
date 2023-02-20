#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <fstream>
#define TINYPLY_IMPLEMENTATION
#include "tinyply.h"

#include "array.h"
#include "hash_table.h"
#include "math_utils.h"
#include "sys_utils.h"

#include "las_normal.h"
#include "las_read.h"
#include "las_source.h"
#include "las_tree.h"

char *get_filename(int x, int y, const char *base_dir)
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
	sprintf(fname + base_len, "%4d_%4d.las", x, y);

	return fname;
}

int inline get_source_idx(const TArray<int> &sources, int source_id)
{
	int n = sources.size;
	while (n--) {
		if (sources[n] == source_id)
			return (n);
	}
	return (-1);
}

inline void rebase_point(struct LasPoint &p, int dx, int dy)
{
	p.x += dx * 100000;
	p.y += dy * 100000;
}

inline bool filter_point(const struct LasPoint &p)
{
	/* 62m of buffer, note that 62.5m = 1/16 km) */
	int bd = 6200;
	bool ret = (p.x > -bd) && (p.x < 100000 + bd) && (p.y > -bd) &&
		   (p.y < 100000 + bd) && (p.classification == 2);
	return (ret);
}

int get_ids_and_counts(TArray<int> &source_ids, TArray<int> &source_counts,
		       int x0, int y0, const char *base_dir)
{
	int file_count = 0;
	for (int i = 0; i < 9; ++i) {
		int dx = (i % 3) - 1;
		int dy = (i / 3) - 1;
		int x = x0 + dx;
		int y = y0 + dy;
		char *fname = get_filename(x, y, base_dir);
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
		 const TArray<int> &counts, TArray<int> &offsets, int x0,
		 int y0, const char *base_dir)
{
	for (int i = 0; i < 9; ++i) {
		int dx = (i % 3) - 1;
		int dy = (i / 3) - 1;
		int x = x0 + dx;
		int y = y0 + dy;
		char *fname = get_filename(x, y, base_dir);
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

void read_and_filter_data(size_t &source_num, TArray<struct LasPoint> &points,
			  int x0, int y0, const char *base_dir)
{
	TArray<int> ids;
	TArray<int> counts;
	get_ids_and_counts(ids, counts, x0, y0, base_dir);
	source_num = ids.size;
	TArray<int> offsets(source_num);
	int point_num = accumulate_counts(offsets, counts);
	points.clear();
	points.resize(point_num);
	fill_points(points, ids, counts, offsets, x0, y0, base_dir);
}

struct CubeTransform {
	float scale;
	float shift_x;
	float shift_y;
	float shift_z;
};

int send_points_to_unit_cube(const TArray<struct LasPoint> &points,
			     struct LasCoords &pos, struct CubeTransform &t)
{
	int min_z = INT_MAX, max_z = -INT_MAX;
	for (size_t i = 0; i < points.size; ++i) {
		min_z = MIN(min_z, points[i].z);
		max_z = MAX(max_z, points[i].z);
	}
	float span = (max_z - min_z) * 0.01f;
	printf("Altitude span : %.0f meters\n", span);
	float n;
	if (span < 1083.33f) {
		n = 6;
	} else if (span < 1625.f) {
		n = 2;
	} else {
		return (-1);
	}
	t.scale = n / (n + 2) * 1e-5;
	t.shift_x = t.shift_y = 1.f / (n + 2);
	float mean = 0.5 * (min_z + max_z) * t.scale;
	t.shift_z = 0.5 - round(16 * mean) * 0.0625;
	assert(min_z * t.scale + t.shift_z >= 0.0625);
	assert(max_z * t.scale + t.shift_z <= 1 - 0.0625);
	for (size_t i = 0; i < points.size; ++i) {
		pos[i].x = points[i].x * t.scale + t.shift_x;
		pos[i].y = points[i].y * t.scale + t.shift_y;
		pos[i].z = points[i].z * t.scale + t.shift_z;
	}
	return (0);
}

int main(int argc, char **argv)
{
	/* Process command line arguments */
	if (argc < 3) {
		printf("Syntax : swiss_lidar x y "
		       "[base_dir]\n");
		return (-1);
	}
	int x0 = atoi(argv[1]);
	int y0 = atoi(argv[2]);
	const char *base_dir = (argc >= 4) ? argv[3] : NULL;
	const int num_probes = (argc >= 5) ? atoi(argv[4]) : 10;

	/* Read data */
	size_t point_num, source_num;
	TArray<struct LasPoint> points;
	read_and_filter_data(source_num, points, x0, y0, base_dir);
	point_num = points.size;
	printf("Total sources : %zu\n", source_num);
	printf("Total filtered points : %zu\n", point_num);

	/* Derive fligh lines azimut from data. Used later for normals
	 */
	printf("Reconstructing flightlines.\n");
	TArray<struct SourceStat> stats(source_num);
	las_stat_sources(points, stats);
	TArray<struct SourceFlightLine> fls(source_num);
	double scale[3] = {1., 1., 1.};
	las_approx_flight_lines(points, scale, stats, fls);

	printf("Rescaling positions to unit cube.\n");
	struct LasCoords pos(point_num);
	struct CubeTransform t;
	if (send_points_to_unit_cube(points, pos, t)) {
		printf("Altitude span too large !\n");
		return (-1);
	}

	/* Build kd-tree for filtered point cloud.
	 */
	printf("Building kd-tree.\n");
	LasTree tree(3, pos, 10);

	/* Build normals using geometry and scanlines*/
	TArray<Vec3> nml(point_num);
	TArray<float> qual(point_num);
	printf("Building unoriented normals :");
	estimate_unoriented_normals(pos, tree, nml, qual, num_probes);
	printf("Orienting normals :\n");
	TArray<EOrient> oriented(point_num, ENone);
	size_t unsettled =
	    orient_normals_with_scanlines(points, fls, qual, nml, oriented);
	printf("    Oriented after scanlines orientation pass : %.0f %%\n",
	       (nml.size - unsettled) * 100.0 / nml.size);
	float progress = 1;
	while (unsettled && progress > 0.01) {
		size_t newly_settled =
		    propagate_normals_once(pos, tree, qual, nml, oriented);
		progress = (float)newly_settled / unsettled;
		unsettled -= newly_settled;
		printf("    Oriented after propagate once pass : %.1f %%\n",
		       (nml.size - unsettled) * 100.0 / nml.size);
	}

	/* Clean-up no longer necessary LAS data */
	points.dispose();
	stats.dispose();
	fls.dispose();

	/* Write oriented points to PLY file */
	size_t final_points = 0;
	for (size_t i = 0; i < point_num; ++i) {
		if (oriented[i] != ENone) {
			pos[final_points] = pos[i];
			nml[final_points] = nml[i];
			qual[final_points] = qual[i];
			final_points++;
		}
	}
	std::filebuf fb_binary;
	char out[16] = {0};
	sprintf(out, "%4d_%4d.ply", x0, y0);
	fb_binary.open(out, std::ios::out | std::ios::binary);
	std::ostream outstream_binary(&fb_binary);
	tinyply::PlyFile file;
	file.add_properties_to_element(
	    "vertex", {"x", "y", "z"}, tinyply::Type::FLOAT32, final_points,
	    reinterpret_cast<uint8_t *>(&pos[0]), tinyply::Type::INVALID, 0);
	file.add_properties_to_element(
	    "vertex", {"nx", "ny", "nz"}, tinyply::Type::FLOAT32, final_points,
	    reinterpret_cast<uint8_t *>(&nml[0]), tinyply::Type::INVALID, 0);
	file.add_properties_to_element(
	    "vertex", {"quality"}, tinyply::Type::FLOAT32, final_points,
	    reinterpret_cast<uint8_t *>(&qual[0]), tinyply::Type::INVALID, 0);
	file.get_comments().push_back("Generated by Didier");
	file.write(outstream_binary, true);

	return (0);
}
