#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>
#include <float.h>

#include <algorithm>
#include <map>
#include <vector>
#include <fstream>

#include "nanoflann.hpp"
#define TINYPLY_IMPLEMENTATION
#include "tinyply.h"

#include "vec3.h"
#include "array.h"
#include "plane_fitting.h"

#include "las_point_cloud.h"
#include "las_source.h"

void print_las_point(const LasPoint &P)
{
	printf("\nPos : %d %d %d\n", P.x, P.y, P.z);
	printf("SourceID : %4d ScanAngle : %+2d GPS : %lf\n", P.source_id,
			P.scan_angle, P.gps_time);
}

struct LasPoint read_las_point(const char *raw, unsigned char point_format)
{
	struct LasPoint P;

	P.x = *(I32 *)(raw + 0);
	P.y = *(I32 *)(raw + 4);
	P.z = *(I32 *)(raw + 8);
	P.intensity = *(U16 *)(raw + 12);
	
	switch (point_format) 
	{
	case 1:
	case 3:
		P.return_number = *(U08 *)(raw + 14) & 0x7;
		P.number_of_returns = *(U08 *)(raw + 14) >> 3 & 0x7;
		P.classification = *(U08 *)(raw + 15);
		P.scan_angle = *(I08 *)(raw + 16);
		P.source_id = *(U16 *)(raw + 18);
		P.gps_time = *(F64 *)(raw + 20);
		break;
	case 6:
	case 7:
	case 8:
		P.return_number = *(U08 *)(raw + 14) & 0x15;
		P.number_of_returns = *(U08 *)(raw + 14) >> 4 & 0x15;
		P.classification = *(U08 *)(raw + 16);
		P.scan_angle = (I08)(*(I16 *)(raw + 18) * 0.006);
		P.source_id = *(U16 *)(raw + 20);
		P.gps_time = *(F64 *)(raw + 22);
		break;

	default:
		assert(0);
	}

	return (P);
}


void print_las_info(const struct LasInfo &info)
{
	printf("LAS Version %d.%d\n", info.version_major, info.version_minor);
	printf("Format: %d, Point Len: %d, Num Points: %zu\n", 
			info.point_format, info.point_size, info.point_num);
	printf("Offsets :\n");
	printf("%lf\n", info.offset[0]);
	printf("%lf\n", info.offset[1]);
	printf("%lf\n", info.offset[2]);
	printf("Scales :\n");
	printf("%lf\n", info.scale[0]);
	printf("%lf\n", info.scale[1]);
	printf("%lf\n", info.scale[2]);
	printf("Bbox :\n");
	char tags[3] = {'x', 'y', 'z'};
	for (int i = 0; i < 3; ++i)
	{
		printf("%c : %lf %lf\n", tags[i], info.min[i], info.max[i]);
	}
}


#define CONSUME(to, bytes) \
	do { \
		if (fread((to), (bytes), 1, f) != 1) return -1; \
		consumed += (bytes); \
	} while (0);

int read_las_info(const char *filename, struct LasInfo &info)
{
	FILE *f;
	f = fopen(filename, "rb");
	if (!f) {
		printf("Las Error: Could not open file %s\n", filename);
		return (-1);
	}

	char buf[64] = {0};
	U16 consumed = 0;

	CONSUME(buf, 4);

	if (strncmp(buf, "LASF", 4) != 0) return -1;

	/* File Source ID */
	CONSUME(buf, 2);
	/* Global Encoding */
	CONSUME(buf, 2);
	/* Project ID 1 */
	CONSUME(buf, 4);
	/* Project ID 2 */
	CONSUME(buf, 2);
	/* Project ID 3 */
	CONSUME(buf, 2);
	/* Project ID 4 */
	CONSUME(buf, 8);
	/* Version Major */
	CONSUME(&info.version_major, 1);
	/* Version Minor */
	CONSUME(&info.version_minor, 1);
	/* System Identifier */
	CONSUME(buf, 32);
	/* Generating Sofware */
	CONSUME(buf, 32);
	/* File Creation Day of Year */
	CONSUME(buf, 2);
	/* File Creation Year */
	CONSUME(buf, 2);

	/* Header size */
	U16 header_size;
	CONSUME(&header_size, 2);
	
	/* Offset to point data */
	CONSUME(&info.offset_to_points, 4);

	/* Number of VLR */
	CONSUME(buf, 4);

	/* Point Data Record Format */
	CONSUME(&info.point_format, 1);
	/* Point Data Record Length */
	CONSUME(&info.point_size, 2);

	U32 legacy_num;
	/* Legacy number of points */
	CONSUME(&legacy_num, 4);
	info.point_num = legacy_num;

	/* Legacy Nbr of Point by Return */
	CONSUME(buf, 20);

	/* Scale */
	CONSUME(info.scale, 24);
	/* Offset */
	CONSUME(info.offset, 24);

	/* Bounding box */
	CONSUME(&info.max[0], 8);
	CONSUME(&info.min[0], 8);
	CONSUME(&info.max[1], 8);
	CONSUME(&info.min[1], 8);
	CONSUME(&info.max[2], 8);
	CONSUME(&info.min[2], 8);

	if (header_size > consumed)
	{
		CONSUME(buf, 8);
		CONSUME(buf, 8);
		CONSUME(buf, 4);
		CONSUME(&info.point_num, 8);
	}

	fclose(f);

	return 0;
}

char *load_raw_las_data(const char *filename, const LasInfo &info,
		        char *dest = NULL)
{
	FILE *f;
	f = fopen(filename, "rb");
	if (!f) return NULL;

	size_t raw_size = info.point_size * info.point_num;

	if (fseek(f, info.offset_to_points, SEEK_SET)) 
	{
		fclose(f);
		return NULL;
	}

	if (!dest)
	{
		dest = (char *)malloc(info.point_size * info.point_num);
		assert(dest);
	}
	
	if (fread(dest, raw_size, 1, f) != 1)
	{
		fclose(f);
		return NULL;
	}

	return dest;
}


int analyze_sources(LasPoint *points, size_t point_num, SourceInfo* sources)
{
	int source_num = 0;
	
	std::map<int, int> id_to_idx;

	for(size_t i = 0; i < point_num; ++i)
	{
		int source_id = points[i].source_id;
		int source_idx;
		auto it = id_to_idx.find(source_id);
		if (it != id_to_idx.end())
		{
			source_idx = (*it).second;
		}
		else
		{
			source_idx = source_num++;
			id_to_idx.insert(std::make_pair(source_id, source_idx));
			sources[source_idx].source_id = source_id;
			sources[source_idx].point_num = 0;
		}
		/* We mutate source ids to source idx */
		points[i].source_id = source_idx;
		sources[source_idx].point_num++;
	}

	printf("Found %d sources.\n", source_num);

	struct SourceStat {
		size_t  idx_of_min = 0;
		size_t  idx_of_max = 0;
		F64  max_gps = -1e9;
		F64  min_gps =  1e9;
		F64  repr_gps;
		F64  tol_gps;
		I08   max_angle = -90;
		I08   min_angle =  90;
		I08   repr_angle;
		bool negative_is_left;
	};

	std::vector<SourceStat> stats(source_num);

	for(size_t i = 0; i < point_num; ++i)
	{
		SourceStat &s = stats[points[i].source_id];
		s.max_angle = std::max(s.max_angle, points[i].scan_angle);
		s.min_angle = std::min(s.min_angle, points[i].scan_angle);
	}

	/* Choose representative scan angles */
	for(int i = 0; i < source_num; ++i)
	{
		SourceStat &s = stats[i];
		int min = s.min_angle;
		int max = s.max_angle;
		if ((max - min) <= 10) 
		{
			s.repr_angle = (max + min) / 2;
		}
		else
		{	
			/* get away by 5° from scan boundaries */
			min += 5;
			max -= 5;
			if (min <= 0 and max >= 0)
			{
				s.repr_angle = 0;
			}
			else if (min >= 0)
			{
				s.repr_angle = min;
			}
			else
			{
				s.repr_angle = max;
			}
		}
		printf("Source %d: min_angle %d max_angle %d repr_angle	%d\n",
				i, s.min_angle, s.max_angle, s.repr_angle);
	}
	
	/* Track start and end of representative flight line */
	bool source_found[source_num];
	bool repr_found[source_num];
	for (int i = 0; i < source_num; ++i)
	{
		source_found[i] = false;
		repr_found[i] = false;
	}
	for(size_t i = 0; i < point_num; ++i)
	{
		int source_id = points[i].source_id;
		if (!source_found[source_id])
		{
			source_found[source_id] = true;
			//printf("Found source %d\n", source_id);
		}
		SourceStat &s = stats[points[i].source_id];
		if (std::abs(points[i].scan_angle - s.repr_angle) > 2)
		{
			continue;
		}
		if (!repr_found[source_id])
		{
			repr_found[source_id] = true;
			//printf("Found representative %d\n", source_id);
		}
		F64 gps_time = points[i].gps_time;
		if (gps_time > s.max_gps)
		{
			s.max_gps = gps_time;
			s.idx_of_max = i;
		}
		if (gps_time < s.min_gps)
		{
			s.min_gps = gps_time;
			s.idx_of_min = i;
		}
	}
	for(int i = 0; i < source_num; ++i)
	{
		//SourceStat &s = stats[i];
		//printf("Times : %lf %lf \n", s.min_gps, s.max_gps);
	}


	/* Set azimuth for sources */
	for(int i = 0; i < source_num; ++i)
	{
		SourceStat &s = stats[i];
		SourceInfo &source = sources[i];
		F64 xe = points[s.idx_of_max].x;
		F64 ye = points[s.idx_of_max].y;
		F64 xs = points[s.idx_of_min].x;
		F64 ys = points[s.idx_of_min].y;
		assert((xe != xs) || (ye != ys));
		source.azimuth = atan2(ye - ys, xe - xs);
	}

	/* Choose representative gps time by source */
	for(int i = 0; i < source_num; ++i)
	{
		SourceStat &s = stats[i];
		s.repr_gps = (s.max_gps + s.min_gps) / 2;
		s.tol_gps = (s.max_gps - s.min_gps) / 100;
		//printf("Times : %lf %lf %lf\n", s.min_gps, s.max_gps, s.tol_gps);
	}

	/* Track left most and right most points at or around repr_gps time */
	for(int i = 0; i < source_num; ++i)
	{
		stats[i].min_angle =  90;
		stats[i].max_angle = -90;
	}
	for(size_t i = 0; i < point_num; ++i)
	{
		SourceStat &s = stats[points[i].source_id];
		F64 gps_time = points[i].gps_time;
		I08 scan_angle = points[i].scan_angle;
		if (std::abs(gps_time - s.repr_gps) > s.tol_gps)
		{
			continue;
		}
		if (scan_angle > s.max_angle)
		{
			s.max_angle  = scan_angle;
			s.idx_of_max = i;
		}
		if (scan_angle < s.min_angle)
		{
			s.min_angle  = scan_angle;
			s.idx_of_min = i;
		}
	}
	
	/* Fill negative_is_left for sources */
	for(int i = 0; i < source_num; ++i)
	{
		SourceStat &s = stats[i];
		SourceInfo &source = sources[i];
		F64 xe = points[s.idx_of_max].x;
		F64 ye = points[s.idx_of_max].y;
		F64 xs = points[s.idx_of_min].x;
		F64 ys = points[s.idx_of_min].y;
		assert((xe != xs) || (ye != ys));
		F64 cross_azimuth = atan2(ye - ys, xe - xs);
		F64 diff_az = source.azimuth - cross_azimuth;
		if (diff_az < 0) diff_az += 2 * M_PI;
		//printf("Source %d : diff_az = %f\n", i, diff_az * 180 / M_PI);
		/* diff_azimuth should be around pi/2 if left is negative
		 * scan_angle and 3pi/2 if right is negative scan_angle */
		// TODO reject source if cross and source mismatch
		source.negative_is_left = diff_az < M_PI;
	}

	return (source_num);
}

struct PointCloudPos {
	size_t point_num;
	const Vec3 *coords;
	inline size_t kdtree_get_point_count() const {return point_num;}
	inline F32  kdtree_get_pt(const size_t idx, int dim) const 
	{
		return coords[idx][dim];
	}
	template <class BBOX>
	bool kdtree_get_bbox(BBOX &bb) const
	{
		(void)bb;
		return false;
	}
};


#define MAX_SOURCES 32
int main(int argc, char **argv)
{
	if (argc < 2) return -1;
	struct LasInfo info;
	
	if (read_las_info(argv[1], info) < 0)
	{
		printf("Error reading %s\n", argv[1]);
	}

	print_las_info(info);

	const char *raw_data = load_raw_las_data(argv[1], info);
	assert(raw_data);

	size_t point_num = info.point_num;

	LasPoint *points = (LasPoint *)malloc(point_num * sizeof(LasPoint));

	for(size_t i = 0; i < point_num; i++)
	{
		const char *raw = raw_data + i * info.point_size;
		points[i] = read_las_point(raw, info.point_format);
	}
	

	struct SourceInfo sources[MAX_SOURCES];
	int source_num = analyze_sources(points, point_num, sources);

	for(int i = 0; i < source_num; ++i)
	{
		SourceInfo &source = sources[i];
		printf("Source %d has azimuth %4d and negative angle is %s\n", 
				source.source_id,
				(int)(source.azimuth * 180 / M_PI),
				source.negative_is_left ? "left" : "right"
				);
	}

	/* Build filter classification remap */
	size_t *remap = (size_t *)malloc(point_num * sizeof(*remap));
	size_t filtered_num = 0;	
	for (size_t i = 0; i < point_num; ++i)
	{
		if (points[i].classification == 2) 
		{
			remap[filtered_num] = i;
			filtered_num++;
		}
	}
	printf("Filtered points : %zu\n", filtered_num);
	if (!filtered_num)
	{
		printf("Exiting.\n");
		exit(0);
	}

	/* Build offseted positions */
	Vec3 *coords  = (Vec3 *)malloc(filtered_num * sizeof(*coords));
	for (size_t i = 0; i < filtered_num; ++i)
	{
		size_t j = remap[i];
		coords[i].x = points[j].x * info.scale[0];
		coords[i].y = points[j].y * info.scale[1];
		coords[i].z = points[j].z * info.scale[2];
	}
	/* Normals */
	Vec3 *normals = (Vec3 *)malloc(filtered_num * sizeof(*normals));

	PointCloudPos cloud = {filtered_num, coords};
	
	using KDTree = nanoflann::KDTreeSingleIndexAdaptor< 
		nanoflann::L2_Simple_Adaptor<float, PointCloudPos>, 
		PointCloudPos, 
		3>;
	
	printf("Building KDTree\n");
	KDTree tree(3, cloud, 10);
	printf("Done.\n");
	
	int num_probes = atoi(argv[2]);
        TArray<unsigned int> ret_index(num_probes);
        TArray<float> out_dist_sqr(num_probes);
	TArray<Vec3>  neigh(num_probes);

	bool *settled = (bool *)malloc(filtered_num * sizeof(*settled));
	F32 *quality = (F32 *)malloc(filtered_num * sizeof(*quality));

	for (size_t i = 0; i < filtered_num; ++i)
	{
		settled[i] = false;
	}

	/* First pass */
	size_t unsettled = 0;
	for (size_t i = 0; i < filtered_num; i++)
	{
		Vec3 query_pt = coords[i];
		int found_probes = tree.knnSearch(&query_pt[0], num_probes, &ret_index[0], 
				&out_dist_sqr[0]);
		if (num_probes != found_probes)
		{
			printf("Discrepancy in probes : %d out of %d\n", found_probes, num_probes);
		}
		for (int j = 0; j < num_probes; ++j)
		{
			neigh[j] = coords[ret_index[j]];
		}
		Vec3 normal;
		//normal = fit_plane(&neigh[0], num_probes);
		F32 qual = fit_plane_vcg<float>(&neigh[0], num_probes, normal);
		quality[i] = qual;
		{
			Vec3 beam;
			F32 scan_angle = points[remap[i]].scan_angle * M_PI / 180;
			SourceInfo &source = sources[points[remap[i]].source_id];
			F32 azimuth = source.azimuth;
			/* For the many .las not following the file spec... */
			if (!source.negative_is_left) scan_angle = -scan_angle;
			beam.x = +sin(azimuth) * sin(scan_angle);
			beam.y = -cos(azimuth) * sin(scan_angle);
			beam.z = -cos(scan_angle);
			F32 test = dot(normal, beam);
			if (std::abs(test) > atof(argv[3]) + 2 * (1 - qual))
			{
				settled[i] = true;
				if (test > 0) normal *= -1;
			}
			else
			{
				unsettled++;
			}
		}
		normals[i] = normal;
		if ((i & 0xFFFFF) == 0) 
		{
			printf("---> %zuM\n", i >> 20);
		}
	}
	printf("Unsettled after init pass : %zu (%f %%)\n", unsettled, 
			(float) unsettled * 100 / filtered_num);

	/* Plague pass */
	size_t newly_settled_num;
	bool *newly_settled = (bool *)malloc(filtered_num * sizeof(*settled));
	do {
		if (argc > 5) break;
		unsettled = 0;
		newly_settled_num = 0;
		for (size_t i = 0; i < filtered_num; i++)
		{
			newly_settled[i] = false;
		}
		for (size_t i = 0; i < filtered_num; i++)
		{
			if (settled[i]) continue;
			Vec3 query_pt = coords[i];
			tree.knnSearch(&query_pt[0], num_probes, &ret_index[0],
					&out_dist_sqr[0]);
			F32 majority = 0;
			for (int j = 0; j < num_probes; ++j)
			{
				int k = ret_index[j];
				if (settled[k]) 
				{
					F32 test = dot(normals[k], normals[i]);
					if (std::abs(test) > (1 + atof(argv[4]) - quality[i]))
					{
						majority += quality[k] * test;
					}
				}
			}
			if (majority != 0) 
			{
				newly_settled[i] = true;
				newly_settled_num++;
				if (majority < 0) normals[i] = -normals[i];
			}
			else
			{
				unsettled++;
			}
		}
		for (size_t i = 0; i < filtered_num; i++)
		{
			settled[i] |= newly_settled[i];
		}
		if (newly_settled_num) 
		{
			printf("Still unsettled : %zu (%f %%)\n", unsettled, 
					(float) unsettled * 100 / filtered_num);
		}
	} while ((float)newly_settled_num / unsettled > 0.001);

	size_t npoints = 0;
	for (size_t i = 0; i < filtered_num; ++i)
	{
		if (settled[i])
		{
			coords[npoints] = coords[i];
			normals[npoints] = normals[i];
			quality[npoints] = quality[i];
			npoints++;
		}
	}

	std::filebuf fb_binary;
	fb_binary.open("normal_test.ply", std::ios::out | std::ios::binary);
	std::ostream outstream_binary(&fb_binary);
	tinyply::PlyFile file;
	file.add_properties_to_element("vertex", { "x", "y", "z" }, 
		tinyply::Type::FLOAT32, npoints, reinterpret_cast<U08*>(coords), tinyply::Type::INVALID, 0);
	file.add_properties_to_element("vertex", { "nx", "ny", "nz" },
		tinyply::Type::FLOAT32, npoints, reinterpret_cast<U08*>(normals), tinyply::Type::INVALID, 0);
	file.add_properties_to_element("vertex", { "quality" },
		tinyply::Type::FLOAT32, npoints, reinterpret_cast<U08*>(quality), tinyply::Type::INVALID, 0);
	file.get_comments().push_back("Generated by Didier");
	file.write(outstream_binary, true);

	free(normals);
	free(coords);
	free(points);
}
