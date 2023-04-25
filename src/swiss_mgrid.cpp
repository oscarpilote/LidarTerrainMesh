#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include "meshoptimizer/src/meshoptimizer.h"

#include "array.h"
#include "hash.h"
#include "hash_table.h"
#include "math_utils.h"
#include "mesh.h"
#include "mesh_ply.h"
#include "mesh_mgrid.h"

struct EncVertex{
	union {
		size_t key;
		struct {
			uint16_t ux;
			uint16_t uy;
			uint16_t uz;
			uint16_t cidx;
		};
	};
};

struct EncVertexHasher {
	static constexpr EncVertex empty = {~static_cast<size_t>(0)};
	size_t hash(EncVertex v) const { return murmur2_64(0, v.key); }
	bool is_empty(EncVertex v) const { return (v.key == empty.key); }
	bool is_equal(EncVertex v1, EncVertex v2) const { return (v1.key == v2.key); }
};

int init_mgrid(const char *fin, const char *fout)
{
	Mesh mesh;
	MBuf data;
	load_ply(mesh, data, fin);

	constexpr uint32_t DIV = 4;
	TArray<CellMeshInfo> info(DIV * DIV);
	TArray<float> shift_x(DIV * DIV);
	TArray<float> shift_y(DIV * DIV);
	TArray<float> max_z(DIV * DIV);
	TArray<float> min_z(DIV * DIV);
	TArray<uint32_t> offset(DIV * DIV);
	/* Initialization */
	for (size_t i = 0; i < DIV * DIV; ++i) {
		info[i].idx_count = 0;
		info[i].vtx_count = 0;
		info[i].scale_z = 2;
		shift_x[i] = (i & (DIV - 1)) * (1.f / DIV);
		shift_y[i] = (i / DIV) * (1.f / DIV);
		max_z[i] = -1;
		min_z[i] = 9;
	}
	/* Split indices accoding to cell of barycenter and set min_z / max_z */
	size_t tri_count = mesh.index_count / 3;
	TArray<uint32_t> tri_idx_to_cell_idx(tri_count);
	const uint32_t *indices = data.indices + mesh.index_offset;
	const Vec3 *positions = data.positions + mesh.vertex_offset;
	for (size_t i = 0; i < tri_count; i++) {
		const Vec3 v0 = positions[indices[3 * i + 0]];
		const Vec3 v1 = positions[indices[3 * i + 1]];
		const Vec3 v2 = positions[indices[3 * i + 2]];
		/* Compute appropriate cell */
		const Vec3 bary = (v0 + v1 + v2) * (1.f / 3.f);
		int cell_x = floor(bary.x * DIV);
		cell_x = (cell_x < 0) ? 0 : cell_x;
		cell_x = (cell_x >= (int)DIV) ? DIV - 1 : cell_x;
		int cell_y = floor(bary.y * DIV);
		cell_y = (cell_y < 0) ? 0 : cell_y;
		cell_y = (cell_y >= (int)DIV) ? DIV - 1 : cell_y;
		int cidx = cell_x + cell_y * DIV;
		/* Update max in min altitudes in cell */
		max_z[cidx] = MAX(max_z[cidx], v0.z);
		max_z[cidx] = MAX(max_z[cidx], v1.z);
		max_z[cidx] = MAX(max_z[cidx], v2.z);
		min_z[cidx] = MIN(min_z[cidx], v0.z);
		min_z[cidx] = MIN(min_z[cidx], v1.z);
		min_z[cidx] = MIN(min_z[cidx], v2.z);
		/* Update dico and idx counts */
		tri_idx_to_cell_idx[i] = cidx;
		info[cidx].idx_count += 3;
	}
	/* Set shift_z */
	for (size_t i = 1; i < DIV * DIV; ++i) {
		printf("Span : %f\n", max_z[i] - min_z[i]);
		float scal = info[i].scale_z;
		info[i].shift_z = scal * floor(100000 * min_z[i] /scal);
	}
	/* Accumulate idx_counts to offset and (tempo.) reset the latter */
	uint32_t total_idx = 0;
	for (size_t i = 0; i < DIV * DIV; ++i) {
		offset[i] = total_idx;
		total_idx += info[i].idx_count;
		info[i].idx_count = 0;
	}

	/* Hash table for vertices in cells */
	HashTable<EncVertex, uint32_t, EncVertexHasher> table(total_idx);
	TArray<uint32_t> cells_idx(total_idx);
	TArray<uint16_t> cells_vtx(total_idx);
	for (size_t i = 0; i < tri_count; ++i) {
		uint32_t cidx = tri_idx_to_cell_idx[i];
		CellMeshInfo &c = info[cidx]; 
		uint32_t v_idx[3];
		for (int k = 0; k < 3; ++k) {
			uint32_t *pval;
			const Vec3 v = positions[indices[3 * i + k]];
			EncVertex ev;
			ev.ux = (v.x - shift_x[cidx]) * (1 << 15) * DIV + (1 << 14);
			ev.uy = (v.y - shift_y[cidx]) * (1 << 15) * DIV + (1 << 14);
			ev.uz = (100000 * v.z - c.shift_z) / c.scale_z;
			ev.cidx = cidx;
			pval = table.get_or_set(ev, c.vtx_count);
			uint16_t *vtx = &cells_vtx[3 * offset[cidx]];
			if (!pval) {
				v_idx[k] = c.vtx_count;
				vtx[c.vtx_count++] = enc - cell_idx;
			} else {
				v_idx[k] = *pval;
			}
		}
		if (v_idx[0] == v_idx[1] || v_idx[1] == v_idx[2] ||
		    v_idx[0] == v_idx[2])
			continue;
		uint32_t *idx = &cells_idx[offset[cidx]];
		idx[c.idx_count++] = v_idx[0];
		idx[c.idx_count++] = v_idx[1];
		idx[c.idx_count++] = v_idx[2];
	}
	/* Encode indices and vertices and write output file */
	FILE *f = fopen(fout, "wb");
	for (size_t i = 0; i < DIV * DIV; ++i) {
		CellInfo &c = cells[i];
		uint32_t *idx = &cells_idx[c.offset];
		meshopt_optimizeVertexCache(idx, idx, c.idx_count, c.vtx_count);
		size_t *vtx = &cells_vtx[c.offset];
		size_t vtx_count = meshopt_optimizeVertexFetch(
		    vtx, idx, c.idx_count, vtx, c.vtx_count, sizeof(size_t));
		printf("%zu %d\n", vtx_count, c.vtx_count);
		// assert(vtx_count == c.vtx_count);
		TArray<uint8_t> vbuf(meshopt_encodeVertexBufferBound(
		    c.vtx_count, sizeof(size_t)));
		vbuf.resize(meshopt_encodeVertexBuffer(
		    &vbuf[0], vbuf.size, vtx, c.vtx_count, sizeof(size_t)));
		TArray<uint8_t> ibuf(
		    meshopt_encodeIndexBufferBound(c.idx_count, c.vtx_count));
		ibuf.resize(meshopt_encodeIndexBuffer(&ibuf[0], ibuf.size, idx,
						      c.idx_count));
		fwrite(vbuf.data, vbuf.size, 1, f);
		fwrite(ibuf.data, ibuf.size, 1, f);
	}
	fclose(f);
	return (0);
}
