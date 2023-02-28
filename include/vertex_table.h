#pragma once

#include <stdint.h>

#include "hash_table.h"
#include "hash.h"
#include "mesh.h"

/* I. Non templated version (uses branching for vertex attribs) */

struct VertexHasher {
	const MBuf* data;
	uint32_t vtx_attr;
	static constexpr uint32_t empty_key = ~static_cast<uint32_t>(0);
	size_t   hash(uint32_t key) const;
	bool is_empty(uint32_t key) const {return (key == empty_key);}
	bool is_equal(uint32_t key1, uint32_t key2) const;
};

using _VertexTable = HashTable<uint32_t, uint32_t, VertexHasher>;

struct VertexTable : public _VertexTable {
	VertexTable(size_t expected_nkeys, const MBuf* data, 
		    uint32_t vtx_attr);
	const MBuf* get_mesh_data() const {return hasher.data;};
	void set_mesh_data(const MBuf* data) {hasher.data = data;}
};

inline
VertexTable::VertexTable(size_t expected_nkeys, const MBuf* data, 
			 uint32_t vtx_attr)
	: _VertexTable::HashTable(expected_nkeys, {data, vtx_attr}) 
{
}

inline size_t
VertexHasher::hash(uint32_t key) const
{
	static_assert(sizeof(uint32_t) == sizeof(float), 
		"uint32_t and float are of different size on this platform.");

	uint32_t hash = 0;
	{
		uint32_t *p = reinterpret_cast<uint32_t*>(data->positions + key);
		hash = murmur2_32(hash, p[0]); 
		hash = murmur2_32(hash, p[1]); 
		hash = murmur2_32(hash, p[2]); 
	}
	if (vtx_attr & VtxAttr::NML)
	{
		uint32_t *p = reinterpret_cast<uint32_t*>(data->normals + key);
		hash = murmur2_32(hash, p[0]); 
		hash = murmur2_32(hash, p[1]); 
		hash = murmur2_32(hash, p[2]); 
	}
	if (vtx_attr & VtxAttr::UV0)
	{
		uint32_t *p = reinterpret_cast<uint32_t*>(data->uv[0] + key);
		hash = murmur2_32(hash, p[0]); 
		hash = murmur2_32(hash, p[1]); 
	}
	if (vtx_attr & VtxAttr::UV1)
	{
		uint32_t *p = reinterpret_cast<uint32_t*>(data->uv[1] + key);
		hash = murmur2_32(hash, p[0]); 
		hash = murmur2_32(hash, p[1]); 
	}
	return (hash);
}

inline bool
VertexHasher::is_equal(uint32_t key1, uint32_t key2) const
{
	bool eq = (data->positions[key1] == data->positions[key2]);

	if (vtx_attr & VtxAttr::NML)
	{
		eq &= (data->normals[key1] == data->normals[key2]); 
	}
	if (vtx_attr & VtxAttr::UV0)
	{
		eq &= (data->uv[0][key1] == data->uv[0][key2]); 
	}
	if (vtx_attr & VtxAttr::UV1)
	{
		eq &= (data->uv[1][key1] == data->uv[1][key2]); 
	}
	return (eq);
}

/* II. Templated version */

template <uint32_t vtx_attr>
struct TVertexHasher {
	const MBuf* data;
	static constexpr uint32_t empty_key = ~static_cast<uint32_t>(0);
	size_t   hash(uint32_t key) const;
	bool is_empty(uint32_t key) const {return (key == empty_key);}
	bool is_equal(uint32_t key1, uint32_t key2) const;
};

template <uint32_t vtx_attr>
using _TVertexTable = HashTable<uint32_t, uint32_t, TVertexHasher<vtx_attr>>;

template <uint32_t vtx_attr>
struct TVertexTable : public _TVertexTable<vtx_attr> {
	TVertexTable(size_t init_num, const MBuf* data); 
};
	
template <uint32_t vtx_attr>
inline 
TVertexTable<vtx_attr>::TVertexTable(size_t expected_nkeys, const MBuf* data)
	: _TVertexTable<vtx_attr>::HashTable(expected_nkeys, {data})
{
}

template <uint32_t vtx_attr>
inline size_t
TVertexHasher<vtx_attr>::hash(uint32_t key) const
{
	static_assert(sizeof(uint32_t) == sizeof(float), 
		"uint32_t and float are of different size on this platform.");

	uint32_t hash = 0;
	{
		uint32_t *p = reinterpret_cast<uint32_t*>(data->positions + key);
		hash = murmur2_32(hash, p[0]); 
		hash = murmur2_32(hash, p[1]); 
		hash = murmur2_32(hash, p[2]); 
	}
	if constexpr ((vtx_attr & VtxAttr::NML) != VtxAttr::None)
	{
		uint32_t *p = reinterpret_cast<uint32_t*>(data->normals + key);
		hash = murmur2_32(hash, p[0]); 
		hash = murmur2_32(hash, p[1]); 
		hash = murmur2_32(hash, p[2]); 
	}
	if constexpr ((vtx_attr & VtxAttr::UV0) != VtxAttr::None)
	{
		uint32_t *p = reinterpret_cast<uint32_t*>(data->uv[0] + key);
		hash = murmur2_32(hash, p[0]); 
		hash = murmur2_32(hash, p[1]); 
	}
	if constexpr ((vtx_attr & VtxAttr::UV1) != VtxAttr::None)
	{
		uint32_t *p = reinterpret_cast<uint32_t*>(data->uv[1] + key);
		hash = murmur2_32(hash, p[0]); 
		hash = murmur2_32(hash, p[1]); 
	}
	return (hash);
}

template <uint32_t vtx_attr>
inline bool
TVertexHasher<vtx_attr>::is_equal(uint32_t key1, uint32_t key2) const
{
	bool eq = (data->positions[key1] == data->positions[key2]);

	if constexpr ((vtx_attr & VtxAttr::NML) != VtxAttr::None)
	{
		eq &= (data->normals[key1] == data->normals[key2]); 
	}

	if constexpr ((vtx_attr & VtxAttr::UV0) != VtxAttr::None)
	{
		eq &= (data->uv[0][key1] == data->uv[0][key2]); 
	}
	
	if constexpr ((vtx_attr & VtxAttr::UV1) != VtxAttr::None)
	{
		eq &= (data->uv[1][key1] == data->uv[1][key2]); 
	}
	
	return (eq);
}
