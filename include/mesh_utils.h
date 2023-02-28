#pragma once

#include <stdint.h>

#include "vec3.h"
#include "aabb.h"
#include "mesh.h"
#include "vertex_table.h"

Aabb compute_mesh_bounds(const Vec3* positions, size_t vertex_count);

Aabb compute_mesh_bounds(const Mesh& mesh, const MBuf& data);

void compute_mesh_normals(const Mesh& mesh, MBuf& data);

void concat_mesh(Mesh& dst_m, MBuf& dst_d, const Mesh& src_m, const MBuf& src_d);

void join_mesh_from_indices(Mesh& dst_m, MBuf& dst_d, const Mesh& src_m, 
		const MBuf& src_d, VertexTable& vtx_table, uint32_t *remap);

void join_mesh_from_vertices(Mesh& dst_m, MBuf& dst_d, const Mesh& src_m, 
		const MBuf& src_d, VertexTable& vtx_table, uint32_t *remap);

void skip_degenerate_tris(Mesh &mesh, MBuf &data);

void compact_mesh(Mesh& mesh, MBuf& data, uint32_t *remap);

void copy_indices(MBuf& dst, size_t dst_off, const MBuf& src, size_t src_off,
		  size_t idx_num, size_t vtx_off = 0);

void copy_vertices(MBuf& dst, size_t dst_off, const MBuf& src, size_t src_off,
		   size_t vtx_num, size_t vtx_off = 0);

uint32_t copy_unique_vertices(MBuf& dst_d, uint32_t dst_off, const MBuf& src_d, 
		uint32_t *vtx_idx, uint32_t vtx_count, VertexTable& vtx_table,
		uint32_t *remap);

