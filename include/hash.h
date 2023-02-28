#pragma once

#include <stdio.h>
#include <stdint.h>

/* Murmur2 32*/
inline uint32_t murmur2_32(uint32_t hash, uint32_t key)
{
	const uint32_t m = 0x5bd1e995;
	const int r = 24;

	key *= m;
	key ^= key >> r;
	key *= m;

	hash *= m;
	hash ^= key;

	return hash;
}

inline uint64_t murmur2_64(uint64_t hash, uint64_t key)
{
	const uint64_t m = 0xc6a4a7935bd1e995llu;
	const int r = 47;

	key *= m;
	key ^= key >> r;
	key *= m;

	hash *= m;
	hash ^= key;

	return hash;
};
	
/* Position hasher */
inline uint32_t position_hash(const float *pos)
{
	const uint32_t *p  = reinterpret_cast<const uint32_t*>(pos);
	unsigned int x = p[0] ^ (p[0] >> 17);
	unsigned int y = p[1] ^ (p[1] >> 17);
	unsigned int z = p[2] ^ (p[2] >> 17);
	return  (x * 73856093u) ^ (y * 19349663u) ^ (z * 83492791u);
}

