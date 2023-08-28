#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "las_point_cloud.h"
#include "las_read.h"

#define CONSUME(to, bytes)                                                     \
	do {                                                                   \
		if (fread((to), (bytes), 1, f) != 1) {                         \
			fclose(f);                                             \
			return -1;                                             \
		}                                                              \
		consumed += (bytes);                                           \
	} while (0);

int las_read_info(const char *filename, struct LasFileInfo &info)
{
	FILE *f = fopen(filename, "rb");
	if (!f)
		return -1;

	char buf[128] = {0};
	uint16_t consumed = 0;

	CONSUME(buf, 4);
	if (strncmp(buf, "LASF", 4) != 0) {
		fclose(f);
		return -1;
	}

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
	int nvlr;
	CONSUME(&nvlr, 4);

	/* Point Data Record Format */
	CONSUME(&info.point_format, 1);
	info.compressed = (info.point_format & 0x80) != 0;
	info.point_format = info.point_format & 0x7F;

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

	if (header_size > consumed) {
		CONSUME(buf, 8);
		CONSUME(buf, 8);
		CONSUME(buf, 4);
		CONSUME(&info.point_num, 8);
		CONSUME(buf, 120);
	}

	assert(header_size == consumed);

	if (info.compressed && header_size == 375 && nvlr >= 1) {
		CONSUME(buf, 54);
		info.copc = strncmp(buf + 2, "copc", 4) == 0;
		info.copc &= *(uint16_t *)(buf + 18) == 1;
	} else {
		info.copc = false;
	}

	fclose(f);
	return 0;
}

void las_print_info(const struct LasFileInfo &info)
{
	printf("LAS Version %d.%d\n", info.version_major, info.version_minor);
	printf("Compressed : %s\n", info.compressed ? "yes" : "no");
	printf("Copc : %s\n", info.copc ? "yes" : "no");
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
	printf("x : %lf %lf\n", info.min[0], info.max[0]);
	printf("y : %lf %lf\n", info.min[1], info.max[1]);
	printf("z : %lf %lf\n", info.min[2], info.max[2]);
}

char *las_load_data(const char *filename, const LasFileInfo &info, char *buf)
{
	FILE *f = fopen(filename, "rb");
	if (!f)
		return NULL;

	if (info.compressed) {
		printf("Cannot read raw data from LAZ file.\n");
		printf("Use copc routines for such files.\n");
		return NULL;
	}

	size_t raw_size = info.point_size * info.point_num;

	if (fseek(f, info.offset_to_points, SEEK_SET)) {
		return NULL;
	}

	bool local_alloc = false;
	if (!buf) {
		buf = (char *)malloc(info.point_size * info.point_num);
		assert(buf);
		local_alloc = true;
	}

	if (fread(buf, raw_size, 1, f) != 1) {
		if (local_alloc)
			free(buf);
		return NULL;
	}

	return buf;
}

struct LasPoint las_read_point(const char *raw, unsigned char point_format)
{
	LasPoint P;

	P.x = *(int32_t *)(raw + 0);
	P.y = *(int32_t *)(raw + 4);
	P.z = *(int32_t *)(raw + 8);
	P.intensity = *(uint16_t *)(raw + 12);

	switch (point_format) {
	case 1:
	case 3:
		P.return_number = *(uint8_t *)(raw + 14) & 0x7;
		P.number_of_returns = *(uint8_t *)(raw + 14) >> 3 & 0x7;
		P.classification = *(uint8_t *)(raw + 15);
		P.scan_angle = *(int8_t *)(raw + 16);
		P.source_id = *(uint16_t *)(raw + 18);
		P.gps_time = *(double *)(raw + 20);
		break;
	case 6:
	case 7:
	case 8:
		P.return_number = *(uint8_t *)(raw + 14) & 0x15;
		P.number_of_returns = *(uint8_t *)(raw + 14) >> 4 & 0x15;
		/* Hack We combine source_id and scanner channel within
		 * source_id (this assumes that source_id does not effectively
		 * uses more than 14bits
		 */
		P.source_id = (*(uint8_t *)(raw + 15) >> 4) << 14;
		P.classification = *(uint8_t *)(raw + 16);
		P.scan_angle = (int8_t)(*(int16_t *)(raw + 18) * 0.006);
		/* Endo of source_id hack */
		P.source_id |= *(uint16_t *)(raw + 20);
		P.gps_time = *(double *)(raw + 22);
		break;

	default:
		assert(0);
	}

	return (P);
}

void print_las_point(const LasPoint &P)
{
	printf("\nPos : %d %d %d\n", P.x, P.y, P.z);
	printf("SourceID : %4d ScanAngle : %+2d GPS : %lf\n", P.source_id,
	       P.scan_angle, P.gps_time);
}
