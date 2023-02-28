#pragma once

#include "mesh.h"

/* I. Non templated version (uses branching for vertex attribs) */
uint32_t build_vertex_remap_old(const Mesh& mesh, const MBuf& data, 
			    uint32_t vtx_attr, uint32_t *remap);

uint32_t build_vertex_remap(const Mesh& mesh, const MBuf& data, 
			    uint32_t vtx_attr, uint32_t *remap);

uint32_t build_vertex_remap_from_indices(const Mesh& mesh, const MBuf& data, 
					 uint32_t vtx_attr, uint32_t *remap);

void remap_index_buffer(const Mesh& mesh, MBuf& data, uint32_t *remap);

/* II. Templated versions */


template <uint32_t vtx_attr>
uint32_t build_vertex_remap(const Mesh& mesh, const MBuf& data, 
			    uint32_t *remap);

template <uint32_t vtx_attr>
uint32_t build_vertex_remap_from_indices(const Mesh& mesh, const MBuf& data, 
					 uint32_t *remap);

constexpr auto build_position_remap = build_vertex_remap<VtxAttr::POS>;

