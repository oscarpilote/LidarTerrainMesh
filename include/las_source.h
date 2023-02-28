#pragma once

#include "array.h"
#include "las_point_cloud.h"

struct SourceStat {
	int point_num;
	double min_gps;
	double max_gps;
	int8_t min_angle;
	int8_t max_angle;
};

struct SourceFlightLine {
	bool is_valid;
	float theta_along;
	float theta_across;
};

int las_get_sources(TArray<struct LasPoint> &points);

void las_stat_sources(const TArray<struct LasPoint> &points,
		      TArray<struct SourceStat> &stats);

int las_approx_flight_lines(const TArray<struct LasPoint> &points,
			    const double *scale,
			    const TArray<struct SourceStat> &stats,
			    TArray<struct SourceFlightLine> &fls);

