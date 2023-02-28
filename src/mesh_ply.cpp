#ifdef DEBUG
#include <stdio.h>
#endif
#include <string.h>

/* For reading fast */
#include "miniply/miniply.h"

/* For writing */
#include <fstream>
#define TINYPLY_IMPLEMENTATION
#include "tinyply/tinyply.h"

#include "mesh.h"
#include "mesh_ply.h"

int load_ply(Mesh &mesh, MBuf &data, const char *fname)
{

	// TODO : should we ?
	data.clear();
	mesh.clear();

	using namespace miniply;
	PLYReader reader(fname);
	if (!reader.valid()) {
		return (EXIT_FAILURE);
	}

	bool got_verts = false;
	bool got_faces = false;

	while (reader.has_element() && (!got_verts || !got_faces)) {
		if (reader.element_is(kPLYVertexElement)) {
			reader.load_element();
			uint32_t vertex_count = reader.num_rows();
			mesh.vertex_offset = 0;
			mesh.vertex_count = vertex_count;

			uint32_t pos_idx[3];
			if (reader.find_pos(pos_idx)) {
				data.vtx_attr |= VtxAttr::POS;
			}

			uint32_t nml_idx[3];
			if (reader.find_normal(nml_idx)) {
				data.vtx_attr |= VtxAttr::NML;
			}

			uint32_t uv_idx[2];
			if (reader.find_texcoord(uv_idx)) {
				data.vtx_attr |= VtxAttr::UV0;
			}

			data.reserve_vertices(vertex_count);

			if (data.vtx_attr & VtxAttr::POS) {
				reader.extract_properties(
				    pos_idx, 3, PLYPropertyType::Float,
				    data.positions);
			}
			if (data.vtx_attr & VtxAttr::NML) {
				reader.extract_properties(
				    nml_idx, 3, PLYPropertyType::Float,
				    data.normals);
			}
			if (data.vtx_attr & VtxAttr::UV0) {
				reader.extract_properties(
				    uv_idx, 2, PLYPropertyType::Float,
				    data.uv[0]);
			}
			got_verts = true;
		} else if (reader.element_is(kPLYFaceElement)) {
			uint32_t idx[1];
			reader.load_element();
			reader.find_indices(idx);
			bool polys = reader.requires_triangulation(idx[0]);
			if (polys && !got_verts) {
				fprintf(stderr,
					"PLY read error in %s: need vertex \
					positions to triangulate faces.\n",
					fname);
				break;
			}
			uint32_t index_count;
			if (polys) {
				index_count = reader.num_triangles(idx[0]) * 3;
				data.reserve_indices(index_count);
				reader.extract_triangles(
				    idx[0], (const float *)data.positions,
				    mesh.vertex_count, PLYPropertyType::Int,
				    data.indices);
			} else {
				index_count = reader.num_rows() * 3;
				data.reserve_indices(index_count);
				reader.extract_list_property(
				    idx[0], PLYPropertyType::Int, data.indices);
			}
			mesh.index_offset = 0;
			mesh.index_count = index_count;
			got_faces = true;
		}
		if (got_verts && got_faces) {
			break;
		}
		reader.next_element();
	}

	if (!got_verts) {
		return (EXIT_FAILURE);
	}

	return (EXIT_SUCCESS);
}

int write_ply(const char *fname, const Mesh &mesh, const MBuf &data)
{
	std::filebuf fbuf;

	fbuf.open(fname, std::ios::out | std::ios::binary);

	std::ostream osb(&fbuf);
	tinyply::PlyFile ply;
	if ((data.vtx_attr & VtxAttr::POS) && mesh.vertex_count) {
		ply.add_properties_to_element(
		    "vertex", {"x", "y", "z"}, tinyply::Type::FLOAT32,
		    mesh.vertex_count,
		    reinterpret_cast<const uint8_t *>(data.positions +
						      mesh.vertex_offset),
		    tinyply::Type::INVALID, 0);
	}
	if ((data.vtx_attr & VtxAttr::NML) && mesh.vertex_count) {
		ply.add_properties_to_element(
		    "vertex", {"nx", "ny", "nz"}, tinyply::Type::FLOAT32,
		    mesh.vertex_count,
		    reinterpret_cast<const uint8_t *>(data.normals +
						      mesh.vertex_offset),
		    tinyply::Type::INVALID, 0);
	}
	if (mesh.index_count) {
		ply.add_properties_to_element(
		    "face", {"vertex_indices"}, tinyply::Type::UINT32,
		    mesh.index_count / 3,
		    reinterpret_cast<const uint8_t *>(data.indices +
						      mesh.index_offset),
		    tinyply::Type::UINT8, 3);
	}

	ply.write(osb, true);

	fbuf.close();
	return (0);
}
