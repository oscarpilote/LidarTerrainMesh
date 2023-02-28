#pragma once

#include <stdint.h>

#include "vec2.h"
#include "vec3.h"

#define MAX_UV_MAPS 2

namespace VtxAttr
{
enum {
	None = 0,
	POS = 1 << 0,
	NML = 1 << 1,
	UV0 = 1 << 2,
	UV1 = 1 << 3,
	MAP = 1 << 4,
	/* Some common combo */
	P = POS,
	PN = POS | NML,
	PNT = POS | NML | UV0,
	PT = POS | UV0,
};
};

/**
 * Mesh Buffer (host side)
 */
struct MBuf {

	uint32_t vtx_attr = 0;

	size_t idx_capacity = 0;
	uint32_t *indices = nullptr;

	size_t vtx_capacity = 0;
	Vec3 *positions = nullptr;
	Vec3 *normals = nullptr;
	Vec2 *uv[MAX_UV_MAPS] = {nullptr};
	uint32_t *remap = nullptr;

	void clear();
	void reserve_indices(size_t num, bool shrink = false);
	void reserve_vertices(size_t num, bool shrink = false);
	void update_vtx_attr(uint32_t attr);
};

/**
 * Mesh (TODO : material)
 */
struct Mesh {
	uint32_t index_offset = 0;
	uint32_t index_count = 0;
	uint32_t vertex_offset = 0;
	uint32_t vertex_count = 0;
	void clear();
};

