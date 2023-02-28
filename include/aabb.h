#pragma once

#include "vec3.h"

template <typename T>
struct TAabb {
	TVec3<T> min;
	TVec3<T> max;

	TAabb<T>& operator|= (const TAabb<T>&& other); 
};
typedef TAabb<float> Aabb;

template<typename T>
TAabb<T>& TAabb<T>::operator|=(const TAabb<T>&& other)
{
	if (other.min.x < min.x) min.x = other.min.x;
	if (other.min.y < min.y) min.y = other.min.y;
	if (other.min.z < min.z) min.z = other.min.z;
	if (other.max.x > max.x) max.x = other.max.x;
	if (other.max.y > max.y) max.y = other.max.y;
	if (other.max.z > max.z) max.z = other.max.z;

	return *this;
}
