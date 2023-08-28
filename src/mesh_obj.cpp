#include <stdio.h>
#include <string.h>

#define FAST_OBJ_IMPLEMENTATION 1
#include "fast_obj/fast_obj.h"
#undef FAST_OBJ_IMPLEMENTATION

#include "array.h"
#include "hash.h"
#include "hash_table.h"
#include "mesh.h"
#include "vec2.h"
#include "vec3.h"

#include "mesh_obj.h"

struct ObjVertexHasher {
	bool has_normals;
	bool has_uv;
	const Vec3 *pos;
	const Vec3 *nml;
	const Vec2 *uv;
	static constexpr fastObjIndex empty_key = {0, 0, 0};
	size_t hash(fastObjIndex key) const;
	bool is_empty(fastObjIndex key) const;
	bool is_equal(fastObjIndex key1, fastObjIndex key2) const;
};
typedef HashTable<fastObjIndex, uint32_t, ObjVertexHasher> ObjVertexTable;

inline size_t ObjVertexHasher::hash(fastObjIndex key) const
{
	return position_hash((const float *)(pos + key.p));
}

bool ObjVertexHasher::is_empty(fastObjIndex key) const { return (key.p == 0); }

bool ObjVertexHasher::is_equal(fastObjIndex key1, fastObjIndex key2) const
{
	bool res = (pos[key1.p] == pos[key2.p]);
	if (res && has_normals)
		res &= (nml[key1.n] == nml[key2.n]);
	if (res && has_uv)
		res &= (uv[key1.t] == uv[key2.t]);

	return (res);
}

static int load_obj(Mesh &mesh, MBuf &data, const fastObjMesh &obj)
{
	data.clear();
	data.vtx_attr = VtxAttr::P;

	size_t obj_vertex_count = 0;
	size_t index_count = 0;

	for (unsigned int i = 0; i < obj.face_count; ++i) {
		obj_vertex_count += obj.face_vertices[i];
		index_count += 3 * (obj.face_vertices[i] - 2);
	}

	data.reserve_indices(index_count);

	bool has_normals = false;
	bool has_uv = false;
	for (size_t i = 0; i < obj_vertex_count; ++i) {
		has_normals |= (obj.indices[i].n != 0);
		has_uv |= (obj.indices[i].t != 0);
	}
	data.vtx_attr |= has_normals ? VtxAttr::NML : 0;
	data.vtx_attr |= has_uv ? VtxAttr::UV0 : 0;

	size_t vertex_count_guess = index_count / 6;
	vertex_count_guess += vertex_count_guess / 2;
	data.reserve_vertices(vertex_count_guess);

	/* Discover vertices and encode indices */

	size_t idx = 0;
	size_t idx_offset = 0;

	size_t vertex_count = 0;
	ObjVertexHasher hasher{has_normals, has_uv, (const Vec3 *)obj.positions,
			       (const Vec3 *)obj.normals,
			       (const Vec2 *)obj.texcoords};
	ObjVertexTable vertices(vertex_count_guess, hasher);

	for (size_t i = 0; i < obj.face_count; ++i) {
		size_t idx_start = idx;
		for (size_t j = 0; j < obj.face_vertices[i]; ++j) {
			if (j >= 3) /* Triangulate */
			{
				data.indices[idx + 0] = data.indices[idx_start];
				data.indices[idx + 1] = data.indices[idx - 1];
				idx += 2;
			}

			fastObjIndex pnt = obj.indices[idx_offset + j];

			uint32_t *p = vertices.get_or_set(pnt, vertex_count);
			if (!p) {
				data.indices[idx] = vertex_count;

				/* Copy vertex */
				if (vertex_count >= data.vtx_capacity) {
					size_t new_cap = data.vtx_capacity +
							 data.vtx_capacity / 2;
					data.reserve_vertices(new_cap);
				}
				void *dst;
				void *src;
				dst = data.positions + vertex_count;
				src = obj.positions + 3 * pnt.p;
				memmove(dst, src, 3 * sizeof(float));
				if (has_normals) {
					dst = data.normals + vertex_count;
					src = obj.normals + 3 * pnt.n;
					memmove(dst, src, 3 * sizeof(float));
				}
				if (has_uv) {
					dst = data.uv[0] + vertex_count;
					src = obj.texcoords + 2 * pnt.t;
					memmove(dst, src, 2 * sizeof(float));
				}
				vertex_count++;
			} else {
				data.indices[idx] = *p;
			}
			idx++;
		}
		idx_offset += obj.face_vertices[i];
	}

	mesh.index_offset = 0;
	mesh.vertex_offset = 0;
	mesh.index_count = index_count;
	mesh.vertex_count = vertex_count;

	bool shrink = true;
	data.reserve_vertices(vertex_count, shrink);

	return (EXIT_SUCCESS);
}

int load_obj(Mesh &mesh, MBuf &data, const char *fname)
{
	fastObjMesh *obj = fast_obj_read(fname);
	if (obj == nullptr) {
		return (EXIT_FAILURE);
	}
	int res = load_obj(mesh, data, *obj);
	fast_obj_destroy(obj);

	return (res);
}

int write_obj(const Mesh &m, const MBuf &data, const char *fname)
{
	FILE *f = fopen(fname, "w");
	if (!f)
		return -1;

	for (size_t i = 0; i < m.vertex_count; i++) {
		const Vec3 &pos = data.positions[i + m.vertex_offset];
		fprintf(f, "v %f %f %f\n", pos.x, pos.y, pos.z);
	}

	if (data.vtx_attr & VtxAttr::NML) {
		for (size_t i = 0; i < m.vertex_count; i++) {
			const Vec3 &nml = data.normals[i + m.vertex_offset];
			fprintf(f, "vn %f %f %f\n", nml.x, nml.y, nml.z);
		}
		for (size_t i = 0; i < m.index_count; i += 3) {
			const uint32_t *idx = data.indices + m.index_offset + i;
			fprintf(f, "f %d//%d %d//%d %d//%d\n", idx[0] + 1,
				idx[0] + 1, idx[1] + 1, idx[1] + 1, idx[2] + 1,
				idx[2] + 1);
		}
	} else {
		for (size_t i = 0; i < m.index_count; i += 3) {
			const uint32_t *idx = data.indices + m.index_offset + i;
			fprintf(f, "f %d %d %d\n", idx[0] + 1, idx[1] + 1,
				idx[2] + 1);
		}
	}
	fclose(f);
	return 0;
}
