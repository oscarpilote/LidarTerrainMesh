#pragma once 

#include "vec3.h"
#include "vec4.h"
#include "quat.h"

/* A ray in 3D space.
 *
 * @param start - starting point of the ray.
 * @param dir - the (positive) direction of the ray.
 */
template <typename T>
struct TRay {
	TVec3<T> start;
	TVec3<T> dir;
};
typedef TRay<float>   Ray;


/* A plane in 3D space.
 *
 * Equation given by
 *
 *      ax + by + cz + d = 0
 */
template <typename T>
struct TPlane {
	union {
		TVec3<T> normal;
		struct {
			T a;
			T b;
			T c;
		};
	};
	T d;
};
typedef TPlane<float> Plane;

/* A sphere in 3D space.
 *
 * @param center - sphere center.
 * @param radius - sphere radius.
 */
template <typename T>
struct TSphere {
	TVec3<T> center;
	T radius;
};
typedef TSphere<float> Sphere;


/* Free functions declarations */

template <typename T>
TPlane<T> plane_from_normal_and_point(const TVec3<T>& normal, const TVec3<T>& point);

template <typename T>
TVec4<T> ray_plane_intersection(const TRay<T>& ray, const TPlane<T>& plane);

template <typename T>
TQuat<T> great_circle_rotation(const TVec3<T>& from, const TVec3<T>& to);

template <typename T>
TVec3<T> normal(const TVec3<T>& v1, const TVec3<T>& v2, const TVec3<T>& v3);

/* Implementations */

template <typename T>
TPlane<T> plane_from_normal_and_point(const TVec3<T>& normal, const TVec3<T>& point)
{
	return {{normal}, {-dot(normal, point)}};
}

template <typename T>
TVec4<T> ray_plane_intersection(const TRay<T>& ray, const TPlane<T>& plane)
{
	T alpha = dot(plane.normal, ray.dir);
	T beta  = -plane.d - dot(plane.normal, ray.start);

	return {alpha * ray.start + beta * ray.dir, alpha};
}	

template <typename T>
TQuat<T> great_circle_rotation(const TVec3<T>& from, const TVec3<T>& to)
{
	assert(approx_equal<T>(norm(from), 1));
	assert(approx_equal<T>(norm(to), 1));

	T cos_angle = dot(from, to);
	if (approx_equal<T>(cos_angle, -1)) 
	{
		return {{0, 0, 1}, 0};
	}

	T cos_half_angle = sqrt((1.f + cos_angle) / 2.f);

	assert(cos_half_angle != 0);
	
	T w = cos_half_angle;
	TVec3<T> xyz = cross(from, to) * (0.5f / cos_half_angle);

	return {xyz, w};
}

template <typename T>
TVec3<T> normal(const TVec3<T>& v1, const TVec3<T>& v2, const TVec3<T>& v3)
{
	return normalized(cross(v2 - v1, v3 - v1));
}

