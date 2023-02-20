#pragma once

#include "las_point_cloud.h"

int las_read_info(const char *filename, struct LasFileInfo &info);

void las_print_info(const struct LasFileInfo &info);

char *las_load_raw_data(const char *filename, const struct LasFileInfo &info,
                        char *buf);

struct LasPoint las_read_point(const char *buf, unsigned char point_format);

void las_print_point(const struct LasPoint &P);
	