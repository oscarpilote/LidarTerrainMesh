#include "mesh.h"

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "sys_utils.h"
#include "vec2.h"
#include "vec3.h"

void MBuf::clear()
{
	MEMFREE(indices);
	idx_capacity = 0;
	MEMFREE(positions);
	MEMFREE(normals);
	MEMFREE(uv[0]);
	MEMFREE(uv[1]);
	MEMFREE(colors);
	MEMFREE(remap);
	vtx_capacity = 0;
}

void MBuf::reserve_indices(size_t num, bool shrink)
{
	assert(num > 0);

	if (num == vtx_capacity)
		return;
	if (num < vtx_capacity && !shrink)
		return;

	REALLOC_NUM(indices, num);
	idx_capacity = num;
}

void MBuf::reserve_vertices(size_t num, bool shrink)
{
	assert(num > 0);

	if (num == vtx_capacity)
		return;
	if (num < vtx_capacity && !shrink)
		return;

	if (vtx_attr & VtxAttr::POS) {
		REALLOC_NUM(positions, num);
	}

	if (vtx_attr & VtxAttr::NML) {
		REALLOC_NUM(normals, num);
	}

	if (vtx_attr & VtxAttr::UV0) {
		REALLOC_NUM(uv[0], num);
	}

	if (vtx_attr & VtxAttr::UV1) {
		REALLOC_NUM(uv[1], num);
	}

	if (vtx_attr & VtxAttr::COL) {
		REALLOC_NUM(colors, num);
	}

	if (vtx_attr & VtxAttr::MAP) {
		REALLOC_NUM(remap, num);
	}

	vtx_capacity = num;
}

void MBuf::add_vtx_attr(uint32_t attr)
{
	uint32_t new_attr = attr & ~vtx_attr;

	vtx_attr |= attr;

	size_t num = vtx_capacity;
	if (!num)
		return;

	if (new_attr & VtxAttr::POS) {
		REALLOC_NUM(positions, num);
	}
	if (new_attr & VtxAttr::NML) {
		REALLOC_NUM(normals, num);
	}
	if (new_attr & VtxAttr::UV0) {
		REALLOC_NUM(uv[0], num);
	}
	if (new_attr & VtxAttr::UV1) {
		REALLOC_NUM(uv[1], num);
	}
	if (new_attr & VtxAttr::COL) {
		REALLOC_NUM(colors, num);
	}
	if (new_attr & VtxAttr::MAP) {
		REALLOC_NUM(remap, num);
	}
}

void MBuf::update_vtx_attr(uint32_t attr)
{
	uint32_t new_attr = attr & ~vtx_attr;
	uint32_t del_attr = vtx_attr & ~attr;

	vtx_attr = attr;

	size_t num = vtx_capacity;
	if (!num)
		return;

	if (new_attr & VtxAttr::POS) {
		REALLOC_NUM(positions, num);
	}
	if (new_attr & VtxAttr::NML) {
		REALLOC_NUM(normals, num);
	}
	if (new_attr & VtxAttr::UV0) {
		REALLOC_NUM(uv[0], num);
	}
	if (new_attr & VtxAttr::UV1) {
		REALLOC_NUM(uv[1], num);
	}
	if (new_attr & VtxAttr::COL) {
		REALLOC_NUM(colors, num);
	}
	if (new_attr & VtxAttr::MAP) {
		REALLOC_NUM(remap, num);
	}

	if (del_attr & VtxAttr::POS) {
		MEMFREE(positions);
	}
	if (del_attr & VtxAttr::NML) {
		MEMFREE(normals);
	}
	if (del_attr & VtxAttr::UV0) {
		MEMFREE(uv[0]);
	}
	if (del_attr & VtxAttr::UV1) {
		MEMFREE(uv[1]);
	}
	if (del_attr & VtxAttr::COL) {
		MEMFREE(colors);
	}
	if (del_attr & VtxAttr::MAP) {
		MEMFREE(remap);
	}
}

void Mesh::clear()
{
	index_offset = 0;
	index_count = 0;
	vertex_offset = 0;
	vertex_count = 0;
}
