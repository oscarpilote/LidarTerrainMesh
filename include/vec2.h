#pragma once

#include <assert.h>
#include <cmath>

template <typename T>
struct TVec2 {
	
	/* Members */
	T x;
	T y;

	/* Constructors */
	constexpr TVec2() = default;
	constexpr TVec2(T x, T y);
	explicit  TVec2(const T* t);

	/* Index Accessor */
	T& operator[] (int n);
	const T& operator[] (int n) const;
	
	/* Equality */
	bool operator== (const TVec2<T>& a) const;

	/* Vector space structure */
	TVec2 operator- () const;
	TVec2& operator+= (const TVec2& a);
	TVec2& operator-= (const TVec2& a);
	TVec2& operator*= (const T& t);
	TVec2& operator/= (const T& t);

	/* Static members */
	static inline TVec2<T> Zero   {0, 0};
	static inline TVec2<T> XAxis  {1, 0};
	static inline TVec2<T> YAxis  {0, 1};
};

typedef TVec2<float> Vec2;

/* Free functions declarations */

template <typename T>
inline TVec2<T> operator+ (const TVec2<T>& a, const TVec2<T>& b);  

template <typename T>
inline TVec2<T> operator- (const TVec2<T>& a, const TVec2<T>& b);  

template <typename T>
inline TVec2<T> operator* (const TVec2<T>& a, const T& t);

template <typename T>
inline TVec2<T> operator* (const T& t, const TVec2<T>& a);

template <typename T>
inline T dot(const TVec2<T>& a, const TVec2<T>& b);

template <typename T>
inline T cross(const TVec2<T>& a, const TVec2<T>& b);

template <typename T>
inline T norm(const TVec2<T> a);

template <typename T>
TVec2<T> normalized(const TVec2<T>& a);

/* Functions implementations */

template <typename T>
inline constexpr TVec2<T>::TVec2(T x, T y): x{x}, y{y} {} 

template <typename T>
inline TVec2<T>::TVec2(const T* t): x{t[0]}, y{t[1]} {}

template <typename T>
inline const T& TVec2<T>::operator[](int n) const
{
	assert(n >= 0 && n <= 1);
	return (&x)[n];
}

template <typename T>
inline T& TVec2<T>::operator[](int n)
{
	assert(n >= 0 && n <= 1);
	return (&x)[n];
}

template <typename T>
inline bool TVec2<T>::operator== (const TVec2<T>& a) const
{
	return (x == a.x && y == a.y);
}

template <typename T>
inline TVec2<T> TVec2<T>::operator-() const
{
	return {-x, -y};
}

template <typename T>
inline TVec2<T>& TVec2<T>::operator+= (const TVec2<T>& a)
{
	x += a.x; 
	y += a.y;
	return (*this);
}

template <typename T>
inline TVec2<T>& TVec2<T>::operator-= (const TVec2<T>& a)
{
	x -= a.x; 
	y -= a.y;
	return (*this);
}

template <typename T>
inline TVec2<T>& TVec2<T>::operator*= (const T& t)
{
	x *= t; 
	y *= t;
	return (*this);
}

template <typename T>
inline TVec2<T>& TVec2<T>::operator/= (const T& t)
{
	x /= t; 
	y /= t;
	return (*this);
}

template <typename T>
inline TVec2<T> operator+ (const TVec2<T>& a, const TVec2<T>& b)
{
	return {a.x + b.x, a.y + b.y};
}

template <typename T>
inline TVec2<T> operator- (const TVec2<T>& a, const TVec2<T>& b)
{
	return {a.x - b.x, a.y - b.y};
}

template <typename T>
inline TVec2<T> operator* (const TVec2<T>& a, const T& t)
{
	return {a.x * t, a.y * t};
}

template <typename T>
inline TVec2<T> operator* (const T& t, const TVec2<T>& a)
{
	return a * t;
}

template <typename T>
inline T dot(const TVec2<T>& a, const TVec2<T>& b)
{
	return (a.x * b.x + a.y * b.y);
}

template <typename T>
inline TVec2<T> cross(const TVec2<T>& a, const TVec2<T>& b)
{
	return (a.x * b.y - a.y * b.x);
}

template <typename T>
inline T norm(const TVec2<T> a)
{
	return sqrt(dot(a,a));
}

template <typename T>
TVec2<T> normalized(const TVec2<T>& a)
{
	T len = norm(a);
	assert(len != 0);
	return a * (1.f / len);
}
