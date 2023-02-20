#pragma once

#include <assert.h>
#include <stdlib.h>

#define ARRAY_FIRST_CAPACITY 8

template <typename T>
struct TArray {
	/**
	 * Members
	 */
	size_t size;
	size_t capacity;
	T *data;
	/**
	 * Methods
	 */
	TArray();
	TArray(const TArray<T> &) = delete;
	TArray(size_t size);
	TArray(size_t size, const T &t);
	~TArray();
	T &operator[](size_t i);
	const T &operator[](size_t i) const;
	void push_back(const T &t);
	T *pop_back();
	void resize(size_t size);
	void reserve(size_t capacity);
	void clear();
	void dispose();
};

template <typename T>
TArray<T>::TArray() : size{0}, capacity{0}, data{nullptr}
{
}

template <typename T>
TArray<T>::TArray(size_t size) : size{size}, capacity{size}
{
	data = static_cast<T *>(malloc(size * sizeof(T)));
}

template <typename T>
TArray<T>::TArray(size_t size, const T &t) : size{size}, capacity{size}
{
	data = static_cast<T *>(malloc(size * sizeof(T)));
	for (size_t i = 0; i < size; i++) {
		data[i] = t;
	}
}

template <typename T>
inline TArray<T>::~TArray()
{
	size = 0;
	capacity = 0;
	free(data);
	data = nullptr;
}

template <typename T>
inline T &TArray<T>::operator[](size_t i)
{
	assert(i < size);
	return (data[i]);
}

template <typename T>
inline const T &TArray<T>::operator[](size_t i) const
{
	assert(i < size);
	return (data[i]);
}

template <typename T>
inline void TArray<T>::push_back(const T &t)
{
	if (size >= capacity) {
		capacity = capacity < ARRAY_FIRST_CAPACITY
			       ? ARRAY_FIRST_CAPACITY
			       : 2 * capacity;
		data = static_cast<T *>(realloc(data, capacity * sizeof(T)));
	}
	data[size++] = t;
}

template <typename T>
inline T *TArray<T>::pop_back()
{
	T *ret = NULL;
	if (size) {
		ret = &data[--size];
	}
	return ret;
}

template <typename T>
void TArray<T>::resize(size_t size)
{
	this->size = size;

	if (size > capacity) {
		data = static_cast<T *>(realloc(data, size * sizeof(T)));
		capacity = size;
	}
}

template <typename T>
void TArray<T>::reserve(size_t capacity)
{
	if (capacity > this->capacity) {
		data = static_cast<T *>(realloc(data, capacity * sizeof(T)));
		this->capacity = capacity;
	}
}

template <typename T>
inline void TArray<T>::clear()
{
	size = 0;
	// capacity = 0;
	free(data);
}

template <typename T>
inline void TArray<T>::dispose()
{
	size = 0;
	capacity = 0;
	free(data);
	data = nullptr;
}
