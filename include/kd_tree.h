#pragma once

#include <stdio.h>

#include "nanoflann/nanoflann.hpp"

#include "array.h"
#include "vec3.h"

/* Cfr nanoflann documentation. */
template <int d>
struct KdCoords {
	const float *pos;
	size_t size;
	inline size_t kdtree_get_point_count() const { return size; }
	inline float kdtree_get_pt(const size_t idx, int dim) const
	{
		assert(idx < size && dim >= 0 && dim < d);
		return pos[d * idx + dim];
	}
	template <class BBOX>
	bool kdtree_get_bbox(BBOX &bb) const
	{
		(void)bb;
		return false;
	}
};

template <int d>
using KdTree = nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<float, KdCoords<d>>, KdCoords<d>, d>;

