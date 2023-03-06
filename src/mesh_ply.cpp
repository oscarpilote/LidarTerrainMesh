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

int load_ply(Mesh &mesh, MBuf &data, const char *fname, int filter)
{

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
			mesh.vertex_count = vertex_count;

			uint32_t extract = 0;
			uint32_t pos_idx[3];
			if ((filter & VtxAttr::POS) &
			    reader.find_pos(pos_idx)) {
				extract |= VtxAttr::POS;
				data.add_vtx_attr(VtxAttr::POS);
			}

			uint32_t nml_idx[3];
			if ((filter & VtxAttr::NML) &
			    reader.find_normal(nml_idx)) {
				extract |= VtxAttr::NML;
				data.add_vtx_attr(VtxAttr::NML);
			}
			uint32_t uv_idx[2];
			if ((filter & VtxAttr::UV0) &
			    reader.find_texcoord(uv_idx)) {
				extract |= VtxAttr::UV0;
				data.add_vtx_attr(VtxAttr::UV0);
			}
			uint32_t col_idx[3];
			if ((filter & VtxAttr::COL) &
			    reader.find_color(col_idx)) {
				extract |= VtxAttr::COL;
				data.add_vtx_attr(VtxAttr::COL);
			}
			data.reserve_vertices(vertex_count +
					      mesh.vertex_offset);

			if (extract & VtxAttr::POS) {
				reader.extract_properties(
				    pos_idx, 3, PLYPropertyType::Float,
				    data.positions + mesh.vertex_offset);
			}
			if (extract & VtxAttr::NML) {
				reader.extract_properties(
				    nml_idx, 3, PLYPropertyType::Float,
				    data.normals + mesh.vertex_offset);
			}
			if (extract & VtxAttr::UV0) {
				reader.extract_properties(
				    uv_idx, 2, PLYPropertyType::Float,
				    data.uv[0] + mesh.vertex_offset);
			}
			if (extract & VtxAttr::COL) {
				reader.extract_properties(
				    col_idx, 3, PLYPropertyType::UChar,
				    data.colors + mesh.vertex_offset);
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
				data.reserve_indices(index_count +
						     mesh.index_offset);
				reader.extract_triangles(
				    idx[0],
				    (const float *)(data.positions +
						    mesh.vertex_offset),
				    mesh.vertex_count, PLYPropertyType::Int,
				    data.indices + mesh.index_offset);
			} else {
				index_count = reader.num_rows() * 3;
				data.reserve_indices(index_count +
						     mesh.index_offset);
				reader.extract_list_property(
				    idx[0], PLYPropertyType::Int,
				    data.indices + mesh.index_offset);
			}
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

int write_ply(const char *fname, const Mesh &mesh, const MBuf &data, int filter)
{
	std::filebuf fbuf;

	fbuf.open(fname, std::ios::out | std::ios::binary);

	std::ostream osb(&fbuf);
	tinyply::PlyFile ply;

	uint32_t extract = filter & data.vtx_attr;

	if ((extract & VtxAttr::POS) && mesh.vertex_count) {
		ply.add_properties_to_element(
		    "vertex", {"x", "y", "z"}, tinyply::Type::FLOAT32,
		    mesh.vertex_count,
		    reinterpret_cast<const uint8_t *>(data.positions +
						      mesh.vertex_offset),
		    tinyply::Type::INVALID, 0);
	}
	if ((extract & VtxAttr::NML) && mesh.vertex_count) {
		ply.add_properties_to_element(
		    "vertex", {"nx", "ny", "nz"}, tinyply::Type::FLOAT32,
		    mesh.vertex_count,
		    reinterpret_cast<const uint8_t *>(data.normals +
						      mesh.vertex_offset),
		    tinyply::Type::INVALID, 0);
	}
	if ((extract & VtxAttr::UV0) && mesh.vertex_count) {
		ply.add_properties_to_element(
		    "vertex", {"u", "v"}, tinyply::Type::FLOAT32,
		    mesh.vertex_count,
		    reinterpret_cast<const uint8_t *>(data.uv[0] +
						      mesh.vertex_offset),
		    tinyply::Type::INVALID, 0);
	}
	if ((extract & VtxAttr::COL) && mesh.vertex_count) {
		ply.add_properties_to_element(
		    "vertex", {"red", "green", "blue"}, tinyply::Type::UINT8,
		    mesh.vertex_count,
		    reinterpret_cast<const uint8_t *>(data.colors +
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
