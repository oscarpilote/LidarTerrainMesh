#pragma once

#include <assert.h>

#include "array.h"
#include "las_tree.h"
#include "vec3.h"

#define NUM_PROBES 10
#define SCAN_TOL 0.25
#define FLIP_TOL 0.35

enum EOrient { ENone = 0, ETmpPlague = 1, EPlague = 2, EScanline = 3 };

void estimate_unoriented_normals(const struct LasCoords &coords,
				 const LasTree &tree, TArray<Vec3> &nml,
				 TArray<float> &qual, int probes);

size_t orient_normals_with_scanlines(const TArray<struct LasPoint> &points,
				     const TArray<struct SourceFlightLine> &fls,
				     TArray<float> &qual, TArray<Vec3> &nml,
				     TArray<EOrient> &oriented,
				     float tol = SCAN_TOL);

size_t propagate_normals_once(const struct LasCoords &pos, const LasTree &tree,
			      TArray<float> &qual, TArray<Vec3> &nml,
			      TArray<EOrient> &oriented,
			      size_t probes = NUM_PROBES, float tol = FLIP_TOL);

