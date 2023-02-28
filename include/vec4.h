#pragma once

#include <assert.h>
#include <cmath>

#include "vec3.h"

template <typename T>
struct alignas(4 * sizeof(T)) TVec4 {
	
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
	constexpr TVec4() = default;
	constexpr TVec4(T x, T y, T z, T w);
	constexpr TVec4(const TVec3<T>& xyz, T w);
	explicit  TVec4(const T* t);
	
	/* Index Accessor */
	T& operator[] (int n);
	const T& operator[] (int n) const;

	/* Equality */
	bool operator== (const TVec4<T>& a) const;


	/* Vector space structure */
	TVec4& operator-  () const;
	TVec4& operator+= (const TVec4& a);
	TVec4& operator-= (const TVec4& a);
	TVec4& operator*= (const T& t);
	TVec4& operator/= (const T& t);

	/* Static members */
	static const TVec4 Zero;
};
typedef TVec4<float> Vec4;

/* Free functions declarations */

template <typename T>
inline TVec4<T> operator+ (const TVec4<T>& a, const TVec4<T>& b);  

template <typename T>
inline TVec4<T> operator- (const TVec4<T>& a, const TVec4<T>& b);  

template <typename T>
inline TVec4<T> operator* (const TVec4<T>& a, const T& t);

template <typename T>
inline TVec4<T> operator* (const T& t, const TVec4<T>& a);

template <typename T>
inline T dot(const TVec4<T>& a, const TVec4<T>& b);

template <typename T>
inline T norm(const TVec4<T> a);

/* Static Helpers */

template <typename T>
const TVec4<T> TVec4<T>::Zero  {0, 0, 0};

/* Functions implementations */

template <typename T>
inline constexpr TVec4<T>::TVec4(T x, T y, T z, T w): x{x}, y{y}, z{z}, w{w} {} 

template <typename T>
inline constexpr TVec4<T>::TVec4(const TVec3<T>& xyz, T w): xyz{xyz}, w{w} {} 

template <typename T>
inline TVec4<T>::TVec4(const T* t): x{t[0]}, y{t[1]}, z{t[2]}, w{t[3]} {}

template <typename T>
inline const T& TVec4<T>::operator[](int n) const
{
	assert(n >= 0 && n <= 3);
	return (&x)[n];
}

template <typename T>
inline T& TVec4<T>::operator[](int n)
{
	assert(n >= 0 && n <= 3);
	return (&x)[n];
}

template <typename T>
inline bool TVec4<T>::operator== (const TVec4<T>& a) const
{
	return (x == a.x && y == a.y && z == a.z && w == a.w);
}

template <typename T>
inline TVec4<T>& TVec4<T>::operator-() const
{
	return {-x, -y, -z, -w};
}

template <typename T>
inline TVec4<T>& TVec4<T>::operator+= (const TVec4<T>& a)
{
	x += a.x; 
	y += a.y;
	z += a.z;
	w += a.w;
	return (*this);
}

template <typename T>
inline TVec4<T>& TVec4<T>::operator-= (const TVec4<T>& a)
{
	x -= a.x; 
	y -= a.y;
	z -= a.z;
	w -= a.w;
	return (*this);
}

template <typename T>
inline TVec4<T>& TVec4<T>::operator*= (const T& t)
{
	x *= t; 
	y *= t;
	z *= t;
	w *= t;
	return (*this);
}

template <typename T>
inline TVec4<T>& TVec4<T>::operator/= (const T& t)
{
	x /= t; 
	y /= t;
	z /= t;
	w /= t;
	return (*this);
}

template <typename T>
inline TVec4<T> operator+ (const TVec4<T>& a, const TVec4<T>& b)
{
	return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
}

template <typename T>
inline TVec4<T> operator- (const TVec4<T>& a, const TVec4<T>& b)
{
	return {a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w};
}

template <typename T>
inline TVec4<T> operator* (const TVec4<T>& a, const T& t)
{
	return {t * a.x, t * a.y, t * a.z, t * a.w};
}

template <typename T>
inline TVec4<T> operator* (const T& t, const TVec4<T>& a)
{
	return a * t;
}

template <typename T>
inline T dot(const TVec4<T>& a, const TVec4<T>& b)
{
	return (a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w);
}

template <typename T>
inline T norm(const TVec4<T> a)
{
	return sqrt(Dot(a,a));
}
