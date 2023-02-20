#pragma once

#include "nanoflann.hpp"

#include "array.h"
#include "vec3.h"

/* Cfr nanoflann documentation. */
struct LasCoords : public TArray<Vec3> {
	LasCoords(size_t size) : TArray<Vec3>(size) {}
	inline size_t kdtree_get_point_count() const { return size; }
	inline float kdtree_get_pt(const size_t idx, int dim) const
	{
		return (*this)[idx][dim];
	}
	template <class BBOX>
	bool kdtree_get_bbox(BBOX &bb) const
	{
		(void)bb;
		return false;
	}
};

using LasTree = nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<float, LasCoords>, LasCoords, 3>;

