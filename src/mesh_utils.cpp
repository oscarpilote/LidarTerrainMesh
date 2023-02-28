#include "mesh_utils.h"

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "aabb.h"
#include "array.h"
#include "geometry.h"
#include "mesh.h"
#include "vec3.h"
#include "vertex_remap.h"
#include "vertex_table.h"

Aabb compute_mesh_bounds(const Vec3 *positions, size_t vertex_count)
{
	Vec3 min = positions[0];
	Vec3 max = positions[0];

	for (size_t i = 1; i < vertex_count; ++i) {
		const Vec3 &pos = positions[i];

		for (size_t j = 0; j < 3; ++j) {
			min[j] = (pos[j] < min[j]) ? pos[j] : min[j];
			max[j] = (pos[j] > max[j]) ? pos[j] : max[j];
		}
	}

	return {min, max};
}

Aabb compute_mesh_bounds(const Mesh &mesh, const MBuf &data)
{
	const Vec3 *positions = data.positions + mesh.vertex_offset;
	size_t vertex_count = mesh.vertex_count;

	return (compute_mesh_bounds(positions, vertex_count));
}

void compute_mesh_normals(const Mesh &mesh, MBuf &data)
{
	if (!(data.vtx_attr & VtxAttr::NML)) {

		void *normals = malloc(data.vtx_capacity * sizeof(Vec3));
		data.normals = static_cast<Vec3 *>(normals);
		data.vtx_attr |= VtxAttr::NML;
	}

	uint32_t *indices = data.indices + mesh.index_offset;
	const Vec3 *positions = data.positions + mesh.vertex_offset;
	Vec3 *normals = data.normals + mesh.vertex_offset;

	/* Init normals to zero */
	for (size_t i = 0; i < mesh.vertex_count; ++i) {
		normals[i] = Vec3::Zero;
	}

	TArray<uint32_t> remap(mesh.vertex_count);
	build_position_remap(mesh, data, &remap[0]);
	// build_vertex_remap_old(mesh, data, VtxAttr::P, &remap[0]);

	for (size_t i = 0; i < mesh.index_count; i += 3) {
		const Vec3 v1 = positions[indices[i + 0]];
		const Vec3 v2 = positions[indices[i + 1]];
		const Vec3 v3 = positions[indices[i + 2]];

		/* Weight normals by triangle area */
		Vec3 n = cross(v2 - v1, v3 - v1);

		/* Accumulate normals of remap targets */
		normals[remap[indices[i + 0]]] += n;
		normals[remap[indices[i + 1]]] += n;
		normals[remap[indices[i + 2]]] += n;
	}

	/* Normalize remap targets and copy them to remap sources */
	for (size_t i = 0; i < mesh.vertex_count; ++i) {
		if (remap[i] == i) {
			normals[i] = normalized(normals[i]);
		} else {
			assert(remap[i] < i);
			normals[i] = normals[remap[i]];
		}
	}
}

void copy_indices(MBuf &dst, size_t dst_off, const MBuf &src, size_t src_off,
		  size_t idx_num, size_t vtx_off)
{
	static_assert(
	    sizeof(*dst.indices) == sizeof(*src.indices),
	    "Error in copy_vertices: dst and src indices type mismatch.");

	assert(src.idx_capacity >= src_off + idx_num);
	assert(dst.idx_capacity >= dst_off + idx_num);

	void *to;
	void *from;

	to = &dst.indices[dst_off];
	from = &src.indices[src_off];
	memcpy(to, from, idx_num * sizeof(*src.indices));

	if (vtx_off) {
		for (size_t i = 0; i < idx_num; ++i) {
			dst.indices[dst_off + i] += vtx_off;
		}
	}
}

void copy_vertices(MBuf &dst, size_t dst_off, const MBuf &src, size_t src_off,
		   size_t vtx_num, size_t vtx_off)
{
	/* Copy only common attributes */
	uint32_t vtx_attr = src.vtx_attr & dst.vtx_attr;

	assert(src.vtx_capacity >= src_off + vtx_num);
	assert(dst.vtx_capacity >= dst_off + vtx_num);

	void *to;
	void *from;

	to = &dst.positions[dst_off];
	from = &src.positions[src_off];
	memcpy(to, from, vtx_num * sizeof(*src.positions));

	if (vtx_attr & VtxAttr::NML) {
		to = &dst.normals[dst_off];
		from = &src.normals[src_off];
		memmove(to, from, vtx_num * sizeof(*src.normals));
	}

	if (vtx_attr & VtxAttr::UV0) {
		to = &dst.uv[0][dst_off];
		from = &src.uv[0][src_off];
		memmove(to, from, vtx_num * sizeof(*src.uv[0]));
	}

	if (vtx_attr & VtxAttr::UV1) {
		to = &dst.uv[1][dst_off];
		from = &src.uv[1][src_off];
		memmove(to, from, vtx_num * sizeof(*src.uv[1]));
	}

	if (vtx_attr & VtxAttr::MAP) {
		to = &dst.remap[dst_off];
		from = &src.remap[src_off];
		memmove(to, from, vtx_num * sizeof(*src.remap));
		if (vtx_off) {
			for (size_t i = 0; i < vtx_num; ++i) {
				dst.remap[dst_off + i] += vtx_off;
			}
		}
	}
}

uint32_t copy_unique_vertices(MBuf &dst_d, uint32_t dst_off, const MBuf &src_d,
			      uint32_t *vtx_idx, uint32_t vtx_count,
			      VertexTable &vtx_table, uint32_t *remap)
{
	/* vtx_table should be based on dst_d and cleared */
	assert(vtx_table.get_mesh_data() == &dst_d && vtx_table.size() == 0);

	/* Source should have all attributes of target */
	assert((dst_d.vtx_attr & src_d.vtx_attr) == dst_d.vtx_attr);

	uint32_t new_vtx_count = 0;
	for (size_t i = 0; i < vtx_count; ++i) {
		size_t vtx_off = dst_off + new_vtx_count;
		copy_vertices(dst_d, vtx_off, src_d, vtx_idx[i], 1);
		uint32_t *p;
		p = vtx_table.get_or_set(vtx_off, new_vtx_count);
		if (p) {
			remap[vtx_idx[i]] = *p;
		} else {
			remap[vtx_idx[i]] = new_vtx_count;
			new_vtx_count++;
		}
	}
	return new_vtx_count;
}

void concat_mesh(Mesh &dst_m, MBuf &dst_d, const Mesh &src_m, const MBuf &src_d)
{
	/* Destination should have no more attributes than source */
	assert((dst_d.vtx_attr & src_d.vtx_attr) == dst_d.vtx_attr);

	size_t total_indices = src_m.index_count + dst_m.index_count;
	dst_d.reserve_indices(dst_m.index_offset + total_indices);

	size_t total_vertices = src_m.vertex_count + dst_m.vertex_count;
	dst_d.reserve_vertices(dst_m.vertex_offset + total_vertices);

	size_t dst_off, src_off, idx_num, vtx_num, vtx_off;

	src_off = src_m.index_offset;
	dst_off = dst_m.index_offset + dst_m.index_count;
	idx_num = src_m.index_count;
	vtx_off = dst_m.vertex_count;
	copy_indices(dst_d, dst_off, src_d, src_off, idx_num, vtx_off);
	dst_m.index_count = total_indices;

	src_off = src_m.vertex_offset;
	dst_off = dst_m.vertex_offset + dst_m.vertex_count;
	vtx_num = src_m.vertex_count;
	vtx_off = dst_m.vertex_count;
	copy_vertices(dst_d, dst_off, src_d, src_off, vtx_num, vtx_off);
	dst_m.vertex_count = total_vertices;
}

void join_mesh_from_indices(Mesh &dst_m, MBuf &dst_d, const Mesh &src_m,
			    const MBuf &src_d, VertexTable &vtx_table,
			    uint32_t *remap)
{
	/* vtx_table should be based on dst_d */
	assert(vtx_table.get_mesh_data() == &dst_d);

	/* Source should have all attributes of target */
	assert((dst_d.vtx_attr & src_d.vtx_attr) == dst_d.vtx_attr);

	uint32_t *idx = dst_d.indices + dst_m.index_offset + dst_m.index_count;
	for (size_t i = 0; i < src_m.index_count; ++i) {
		size_t dst_off = dst_m.vertex_offset + dst_m.vertex_count;
		size_t src_idx = src_d.indices[src_m.index_offset + i];
		size_t src_off = src_idx + src_m.vertex_offset;
		copy_vertices(dst_d, dst_off, src_d, src_off, 1);
		uint32_t *p;
		p = vtx_table.get_or_set(dst_off, dst_m.vertex_count);
		if (p) {
			idx[i] = *p;
		} else {
			idx[i] = dst_m.vertex_count;
			dst_m.vertex_count++;
		}
		if (remap) {
			assert(src_idx < src_m.vertex_count);
			remap[src_idx] = idx[i];
		}
	}
	dst_m.index_count += src_m.index_count;
}

void join_mesh_from_vertices(Mesh &dst_m, MBuf &dst_d, const Mesh &src_m,
			     const MBuf &src_d, VertexTable &vtx_table,
			     uint32_t *remap)
{
	/* vtx_table should be based on dst_d */
	assert(vtx_table.get_mesh_data() == &dst_d);

	/* Source should have all attributes of target */
	assert((dst_d.vtx_attr & src_d.vtx_attr) == dst_d.vtx_attr);

	for (size_t i = 0; i < src_m.vertex_count; ++i) {
		size_t dst_off = dst_m.vertex_offset + dst_m.vertex_count;
		size_t src_off = src_m.vertex_offset + i;
		copy_vertices(dst_d, dst_off, src_d, src_off, 1);
		uint32_t *p;
		p = vtx_table.get_or_set(dst_off, dst_m.vertex_count);
		if (p) {
			remap[i] = *p;
		} else {
			remap[i] = dst_m.vertex_count;
			dst_m.vertex_count++;
		}
	}

	uint32_t *dst_idx =
	    dst_d.indices + dst_m.index_offset + dst_m.index_count;
	uint32_t *src_idx = src_d.indices + src_m.index_offset;
	for (size_t j = 0; j < src_m.index_count; ++j) {
		dst_idx[j] = remap[src_idx[j]];
	}
	dst_m.index_count += src_m.index_count;
}

void compact_mesh(Mesh &mesh, MBuf &data, uint32_t *remap)
{
	/* Eliminate spatially degenerate triangles */
	build_position_remap(mesh, data, remap);
	uint32_t total_indices = 0;
	uint32_t *idx = data.indices + mesh.index_offset;
	for (size_t i = 0; i < mesh.index_count; i += 3) {
		uint32_t i1 = remap[idx[i + 0]];
		uint32_t i2 = remap[idx[i + 1]];
		uint32_t i3 = remap[idx[i + 2]];
		if (i1 == i2 || i1 == i3 || i2 == i3) {
			continue;
		}
		idx[total_indices + 0] = idx[i + 0];
		idx[total_indices + 1] = idx[i + 1];
		idx[total_indices + 2] = idx[i + 2];
		total_indices += 3;
	}
	mesh.index_count = total_indices;

	// remap_index_buffer(mesh, data, remap);
	// remap_vertex_buffer(mesh, data, remap);
	// mesh.vertex_count = vtx_num;
}

void skip_degenerate_tris(Mesh &mesh, MBuf &data)
{
	uint32_t new_idx_count = 0;
	uint32_t *idx = data.indices + mesh.index_offset;
	for (uint32_t k = 0; k < mesh.index_count; k += 3) {
		uint32_t i0 = idx[k + 0];
		uint32_t i1 = idx[k + 1];
		uint32_t i2 = idx[k + 2];
		if (i0 == i1 || i1 == i2 || i0 == i2)
			continue;
		idx[new_idx_count++] = i0;
		idx[new_idx_count++] = i1;
		idx[new_idx_count++] = i2;
	}
	mesh.index_count = new_idx_count;
}
