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



struct LasInfo  {
	uint64_t point_num;
	double   offset[3];
	double   scale[3];
	double   min[3];
	double   max[3];
	uint8_t  version_major;
	uint8_t  version_minor;
	uint8_t  point_format;
	uint32_t offset_to_points;
	uint16_t point_size;
};


struct LasPoint {
	int32_t  x;
	int32_t  y;
	int32_t  z;
	uint16_t intensity;
	uint8_t  return_number;
	uint8_t  number_of_returns;
	uint16_t classification;
	 int8_t  scan_angle;
	 int8_t  unused;
	uint16_t source_id;
	double   gps_time;
};

void print_las_point(LasPoint &P)
{
	printf("\nPos : %d %d %d\n", P.x, P.y, P.z);
	printf("SourceID : %4d ScanAngle : %+2d GPS : %lf\n", P.source_id,
			P.scan_angle, P.gps_time);
}


LasPoint read_point(const char *raw, unsigned char point_format)
{
	LasPoint P;

	P.x = *(int32_t *)(raw + 0);
	P.y = *(int32_t *)(raw + 4);
	P.z = *(int32_t *)(raw + 8);
	P.intensity = *(uint16_t *)(raw + 12);
	switch (point_format) 
	{
	case 1:
		P.return_number = *(uint8_t *)(raw + 14) & 0x7;
		P.number_of_returns = *(uint8_t *)(raw + 14) >> 3 & 0x7;
		P.classification = *(uint8_t *)(raw + 15);
		P.scan_angle = *(int8_t *)(raw + 16);
		P.source_id = *(uint16_t *)(raw + 18);
		P.gps_time = *(double *)(raw + 20);
		break;
	case 6:
		P.return_number = *(uint8_t *)(raw + 14) & 0x15;
		P.number_of_returns = *(uint8_t *)(raw + 14) >> 4 & 0x15;
		P.classification = *(uint8_t *)(raw + 16);
		P.scan_angle = (int8_t)(*(int16_t *)(raw + 18) * 0.006);
		P.source_id = *(uint16_t *)(raw + 20);
		P.gps_time = *(double *)(raw + 22);
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
	uint16_t consumed = 0;

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
	uint16_t header_size;
	CONSUME(&header_size, 2);
	
	/* Offset to point data */
	CONSUME(&info.offset_to_points, 4);

	/* Number of VLR */
	CONSUME(buf, 4);

	/* Point Data Record Format */
	CONSUME(&info.point_format, 1);
	/* Point Data Record Length */
	CONSUME(&info.point_size, 2);

	uint32_t legacy_num;
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

struct PointCloudPos {
	size_t point_num;
	const Vec3 *coords;
	inline size_t kdtree_get_point_count() const {return point_num;}
	inline float  kdtree_get_pt(const size_t idx, int dim) const 
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
		points[i] = read_point(raw, info.point_format);
	}

	free((void *)raw_data);

	printf("Got LIDAR.\n");
	
	/* Build filter classification remap */
	size_t *remap = NULL;
	if (argc > 3)
	{
		int classif = atoi(argv[3]);
		remap = (size_t *)malloc(point_num * sizeof(*remap));
		size_t filtered_num = 0;	
		for (size_t i = 0; i < point_num; ++i)
		{
			if (points[i].classification == classif) 
			{
				remap[filtered_num] = i;
				filtered_num++;
			}
		}
		point_num = filtered_num;
	}


	/* Build offseted positions */
	Vec3 *coords  = (Vec3 *)malloc(point_num * sizeof(*coords));
	for (size_t i = 0; i < point_num; ++i)
	{
		size_t j = remap ? remap[i] : i;
		coords[i].x = points[j].x * info.scale[0];
		coords[i].y = points[j].y * info.scale[1];
		coords[i].z = points[j].z * info.scale[2];
	}
	/* Normals */
	Vec3 *normals = (Vec3 *)malloc(point_num * sizeof(*normals));

	PointCloudPos cloud = {point_num, coords};
	
	using KDTree = nanoflann::KDTreeSingleIndexAdaptor< 
		nanoflann::L2_Simple_Adaptor<float, PointCloudPos>, 
		PointCloudPos, 
		3>;
	
	printf("Building KDTree.\n");
	KDTree tree(3, cloud, 10);
	printf("KDTree ready.\n");
	
	//int num_probes = atoi(argv[2]);
	float radius = atof(argv[2]);
        //TArray<unsigned int> ret_index(num_probes);
        //TArray<float> out_dist_sqr(num_probes);
	//TArray<Vec3>  neigh(num_probes);
	
	std::vector<Vec3> neigh;
	std::vector<std::pair<uint32_t, float>> ret_matches;
        nanoflann::SearchParams params;
        params.sorted = false;

	float *qual1 = (float *)malloc(point_num * sizeof(*qual1));
	float *qual2 = (float *)malloc(point_num * sizeof(*qual2));

	/* Quality test */
	for (size_t i = 0; i < point_num; i++)
	{
		Vec3 query_pt = coords[i];
		const size_t neigh_num = tree.radiusSearch(
			&query_pt[0], radius, ret_matches, params);
		neigh.resize(neigh_num);
		for (size_t j = 0; j < neigh_num; ++j)
		{
			neigh[j] = coords[ret_matches[j].first];
		}
		fit_plane_qual<float>(&neigh[0], neigh_num, normals[i],
				qual1[i], qual2[i]);
		if (normals[i].z < 0) normals[i] = -normals[i];
		if ((i & 0xFFFFF) == 0) 
		{
			printf("---> %zuM\n", i >> 20);
		}
	}

	std::filebuf fb_binary;
	fb_binary.open("quality_test.ply", std::ios::out | std::ios::binary);
	std::ostream outstream_binary(&fb_binary);
	tinyply::PlyFile file;
	file.add_properties_to_element("vertex", { "x", "y", "z" }, 
		tinyply::Type::FLOAT32, point_num,
		reinterpret_cast<uint8_t*>(coords), tinyply::Type::INVALID, 0);
	file.add_properties_to_element("vertex", { "nx", "ny", "nz" },
		tinyply::Type::FLOAT32, point_num,
		reinterpret_cast<uint8_t*>(normals), tinyply::Type::INVALID, 0);
	file.add_properties_to_element("vertex", { "qual1" },
		tinyply::Type::FLOAT32, point_num, 
		reinterpret_cast<uint8_t*>(qual1), tinyply::Type::INVALID, 0);
	file.add_properties_to_element("vertex", { "qual2" },
		tinyply::Type::FLOAT32, point_num, 
		reinterpret_cast<uint8_t*>(qual2), tinyply::Type::INVALID, 0);
	file.get_comments().push_back("Generated by Didier");
	file.write(outstream_binary, true);

	free(remap);
	free(normals);
	free(coords);
	free(points);
}
