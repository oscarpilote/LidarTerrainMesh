#include "vertex_remap.h"

#include <stdint.h>
#include <stdio.h>

#include "mesh.h"
#include "vertex_table.h"


/* I. Non templated version (uses branching for vertex attribs) */

uint32_t build_vertex_remap_old(const Mesh& mesh, const MBuf& data, 
			    uint32_t vtx_attr, uint32_t *remap)
{
	VertexTable vtx_table{mesh.vertex_count, &data, vtx_attr};

	uint32_t num = 0;

	for (size_t i = 0; i < mesh.vertex_count; ++i)
	{
		uint32_t *p = vtx_table.get_or_set(i + mesh.vertex_offset, i);
		
		if (p)
		{	
			remap[i] = *p;
		}
		else
		{
			remap[i] = i;
			num++;
		}
	}
	return (num);
}

uint32_t build_vertex_remap(const Mesh& mesh, const MBuf& data, 
			    uint32_t vtx_attr, uint32_t *remap)
{

	switch (vtx_attr) {
	case (VtxAttr::P):
		return build_vertex_remap<VtxAttr::P>(mesh, data, remap);
	case (VtxAttr::PN):
		return build_vertex_remap<VtxAttr::PN>(mesh, data, remap);
	case (VtxAttr::PNT):
		return build_vertex_remap<VtxAttr::PNT>(mesh, data, remap);
	case (VtxAttr::PT):
		return build_vertex_remap<VtxAttr::PT>(mesh, data, remap);
	default:
		assert(0);
		return 0;
	}
}

uint32_t build_vertex_remap_from_indices(const Mesh& mesh, const MBuf& data, 
					 uint32_t vtx_attr, uint32_t *remap)
{
	VertexTable vtx_table{mesh.vertex_count, &data, vtx_attr};

	/* Init remap with ~0 == NO-REMAP */
	for (size_t i = 0; i < mesh.vertex_count; ++i)
	{
		remap[i] = ~static_cast<uint32_t>(0);
	}

	uint32_t num = 0;
	for (size_t i = 0; i < mesh.index_count; ++i)
	{
		uint32_t idx = data.indices[mesh.index_offset + i];
		uint32_t key = mesh.vertex_offset + idx;
		uint32_t *p = vtx_table.get_or_set(key, num);
		
		if (p)
		{
			remap[idx] = *p;
		}
		else
		{
			remap[idx] = num;
			num++;
		}
	}
	return (num);
}


/* Templated version */


template <uint32_t vtx_attr>
uint32_t build_vertex_remap(const Mesh& mesh, const MBuf& data, uint32_t *remap)
{
	TVertexTable<vtx_attr> vtx_remap {mesh.vertex_count, {&data}};

	uint32_t num = 0;

	for (size_t i = 0; i < mesh.vertex_count; ++i)
	{
		uint32_t *p = vtx_remap.get_or_set(i + mesh.vertex_offset, i);
		
		if (p)
		{
			remap[i] = *p;
		}
		else
		{
			remap[i] = i;
			num++;
		}
	}
	return (num);
}

template <uint32_t vtx_attr>
uint32_t build_vertex_remap_from_indices(const Mesh& mesh, const MBuf& data, 
					 uint32_t *remap)
{
	TVertexTable<vtx_attr> vtx_remap {mesh.vertex_count, {&data}};

	/* Init remap with ~0 == NO-REMAP */
	for (size_t i = 0; i < mesh.vertex_count; ++i)
	{
		remap[i] = ~static_cast<uint32_t>(0);
	}

	uint32_t num = 0;
	for (size_t i = 0; i < mesh.index_count; ++i)
	{
		uint32_t idx = data.indices[mesh.index_offset + i];
		uint32_t key = mesh.vertex_offset + idx;
		uint32_t *p = vtx_remap.get_or_set(key, num);
		
		if (p)
		{
			remap[idx] = *p;
		}
		else
		{
			remap[idx] = num;
			num++;
		}
	}
	return (num);
}

/* Instantiations */

template uint32_t build_vertex_remap<VtxAttr::P>(const Mesh& mesh, 
		const MBuf& data, uint32_t *remap);

template uint32_t build_vertex_remap_from_indices<VtxAttr::P>(
		const Mesh& mesh, const MBuf& data, uint32_t *remap);

template uint32_t build_vertex_remap<VtxAttr::PN>(
		const Mesh& mesh, const MBuf& data, uint32_t *remap);

template uint32_t build_vertex_remap_from_indices<VtxAttr::PN>(
		const Mesh& mesh, const MBuf& data, uint32_t *remap);

template uint32_t build_vertex_remap<VtxAttr::PT>(
		const Mesh& mesh, const MBuf& data, uint32_t *remap);

template uint32_t build_vertex_remap_from_indices<VtxAttr::PT>(
		const Mesh& mesh, const MBuf& data, uint32_t *remap);

template uint32_t build_vertex_remap<VtxAttr::PNT>(
		const Mesh& mesh, const MBuf& data, uint32_t *remap);

template uint32_t build_vertex_remap_from_indices<VtxAttr::PNT>(
		const Mesh& mesh, const MBuf& data,	uint32_t *remap);


/* Apply remaps */

void remap_index_buffer(const Mesh& mesh, MBuf& data, uint32_t *remap)
{
	uint32_t *idx = data.indices + mesh.index_offset;
	for (size_t i = 0; i < mesh.index_count; ++i)
	{
		assert(idx[i] < mesh.vertex_count);
		idx[i] = remap[idx[i]];
	}
}


