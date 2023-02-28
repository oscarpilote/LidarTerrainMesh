#pragma once

#include <assert.h>

#include "math_utils.h"
#include "vec3.h"

template <typename T>
struct alignas(4 * sizeof(T)) TQuat {
	
	/* Members */
	union {
		TVec3<T> xyz;
		struct {
			T x;
			T y;
			T z;
		};
	};
	T w;

	/* Constructors */
	TQuat() = default;
	constexpr TQuat(T x, T y, T z, T w);
	constexpr TQuat(TVec3<T> xyz, T w);

	/* Static Helpers */
	static inline TQuat Identity {TVec3<T>::Zero, 1};

	/* Operator and Methods */
	TQuat& operator*= (const TQuat& a);
	TQuat& operator*= (const T& t);
	TQuat conj() const;
	TQuat inv() const;
	TQuat& normalise();
};

typedef TQuat<float> Quat;

/* Free functions */

/* Note: -q := q.conj() , this is handy for when working with
 * unit quaternions as rotations because then -q is the unit
 * quaternion associated to the inverse rotation. 
 * In brief: if q is rotate then -q is unrotate.
 */
template <typename T>
inline TQuat<T> operator- (const TQuat<T>& a);

template <typename T>
inline TQuat<T> operator* (const TQuat<T>& a, const TQuat<T>& b);

template <typename T>
inline TQuat<T> operator* (const TQuat<T>& a, const T& t);

template <typename T>
inline TQuat<T> operator* (const T& t, const TQuat<T>& a);

template <typename T>
T dot(const TQuat<T>& a, const TQuat<T>& b);

template <typename T>
inline T norm(const TQuat<T>& a);

template <typename T>
TQuat<T> pow(TQuat<T>& q, T t);

template <typename T>
TQuat<T> slerp(TQuat<T>& q0, TQuat<T>q1, T t);

/* Methods implementation */

template <typename T>
inline constexpr TQuat<T>::TQuat(T x, T y, T z, T w): x{x}, y{y}, z{z}, w{w} {}

template <typename T>
inline constexpr TQuat<T>::TQuat(TVec3<T> xyz, T w): xyz{xyz}, w{w} {}

template <typename T>
inline TQuat<T>& TQuat<T>::operator*= (const TQuat& a)
{
	/* TODO : better ? */
	*this = *this * a;
	return (*this);
}

template <typename T>
inline TQuat<T>& TQuat<T>::operator*= (const T& t)
{
	xyz *= t;
	w *= t;
	return (*this);
}

template <typename T>
inline TQuat<T> TQuat<T>::conj() const
{
	return {-xyz, w};
}

template <typename T>
inline TQuat<T> TQuat<T>::inv() const
{
	return this->conj() *= (1.f / dot(*this, *this));
}

template <typename T>
inline TQuat<T>& TQuat<T>::normalise()
{
	*this *= (1.f / norm(*this));
	return (*this);
}

/* Free functions implementation */

template <typename T>
inline TQuat<T> operator- (const TQuat<T>& a)
{
	return a.conj();
}

template <typename T>
inline TQuat<T> operator* (const TQuat<T>& a, const TQuat<T>& b)
{
	return {a.w * b.xyz + a.xyz * b.w + cross(a.xyz, b.xyz),
		a.w * b.w - dot(a.xyz, b.xyz)};
}

template <typename T>
inline TQuat<T> operator* (const TQuat<T>& a, const T& t)
{
	return {t * a.xyz, t * a.w};
}

template <typename T>
inline TQuat<T> operator* (const T& t, const TQuat<T>& a)
{
	return a * t;
}

template <typename T>
inline T dot(const TQuat<T>& a, const TQuat<T>& b)
{
	return dot(a.xyz, b.xyz) + a.w * b.w;
}

template <typename T>
inline T norm(const TQuat<T>& a)
{
	return sqrt(dot(a, a));
}

template <typename T>
TQuat<T> pow(TQuat<T>& q, T t)
{
	T qw = q.w < -1 ? -1 : q.w > 1 ? 1 : q.w;
	T omega = acos(qw);
	T cto = cos(t * omega);
	T sto = sin(t * omega);
	T so  = sin(omega);
	T mult = (so != 0) ? sto / so : t;

	return {mult * q.xyz, cto};
}

template <typename T>
TQuat<T> slerp(TQuat<T>& q0, TQuat<T>q1, T t)
{
	return (q0 * pow(q0.inv() * q1, t));
}




