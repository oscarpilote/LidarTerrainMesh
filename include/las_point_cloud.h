#pragma once

#include <stdint.h>

struct LasFileInfo
{
	uint64_t point_num;
	double offset[3];
	double scale[3];
	double min[3];
	double max[3];
	uint8_t version_major;
	uint8_t version_minor;
	uint8_t point_format;
	uint32_t offset_to_points;
	uint16_t point_size;
};

struct LasPoint
{
	int32_t x;		   //  0 -  4
	int32_t y;		   //  4 -  8
	int32_t z;		   //  8 - 12
	uint16_t intensity;	   // 12 - 14
	uint8_t return_number;	   // 14 - 15
	uint8_t number_of_returns; // 15 - 16
	uint16_t classification;   // 16 - 18
	int8_t scan_angle;	   // 18 - 19
	int8_t pad1;		   // 19 - 20
	uint16_t source_id;	   // 20 - 22
	uint16_t source_idx;	   // 22 - 24
	double gps_time;	   // 24 - 32
};
