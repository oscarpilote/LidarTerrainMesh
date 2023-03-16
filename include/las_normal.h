#pragma once

#include "kd_tree.h"
#include "vec3.h"

#include "las_source.h"

#define ORIENT_PROBES 10
#define SCAN_TOL 0.25
#define FLIP_PROBES 6
#define FLIP_TOL 0.85
#define NML_Z_THRESH 0.55

enum EOrient {
	ENone = 0,
	ETmpPlague = 1,
	EPlague = 2,
	EPositiveZ = 3,
	EScanline = 4
};

void estim_unoriented_nml(const Vec3 *pos, size_t point_num, Vec3 *nml,
			  float *qual, const KdTree<3> &tree,
			  void (*cb)(float progress) = NULL,
			  int probes = ORIENT_PROBES);

size_t orient_nml_with_z(Vec3 *nml, EOrient *oriented, size_t point_num,
			 const float *qual, float tol = NML_Z_THRESH);

size_t orient_nml_with_scan(const LasPoint *points, size_t point_num,
			    const SourceFlightLine *fls, size_t source_num,
			    const float *qual, Vec3 *nml, EOrient *oriented,
			    float tol = SCAN_TOL);

size_t propagate_nml_once(const Vec3 *pos, size_t point_num,
			  const KdTree<3> &tree, const float *qual, Vec3 *nml,
			  EOrient *oriented, size_t probes = FLIP_PROBES,
			  float tol = FLIP_TOL);
