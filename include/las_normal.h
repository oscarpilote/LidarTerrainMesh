#pragma once

#include "vec3.h"

#include "las_source.h"
#include "las_tree.h"

#define NUM_PROBES 10
#define SCAN_TOL 0.25
#define FLIP_TOL 0.35
#define NML_Z_THRESH 0.3

enum EOrient {
	ENone = 0,
	ETmpPlague = 1,
	EPlague = 2,
	EScanline = 3,
	EPositiveZ = 4
};

void estim_unoriented_nml(const Vec3 *pos, size_t point_num, Vec3 *nml,
			  float *qual, const LasTree &tree,
			  void (*cb)(float progress) = NULL,
			  int probes = NUM_PROBES);

size_t orient_nml_with_z(Vec3 *nml, EOrient *oriented, size_t point_num);

size_t orient_nml_with_scan(const LasPoint *points, size_t point_num,
			    const SourceFlightLine *fls, size_t source_num,
			    const float *qual, Vec3 *nml, EOrient *oriented,
			    float tol = SCAN_TOL);

size_t propagate_nml_once(const Vec3 *pos, size_t point_num,
			  const LasTree &tree, const float *qual, Vec3 *nml,
			  EOrient *oriented, size_t probes = NUM_PROBES,
			  float tol = FLIP_TOL);