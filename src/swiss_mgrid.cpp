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

struct EncVertexHasher {
	static constexpr size_t empty_key = ~static_cast<size_t>(0);
	size_t hash(size_t key) const { return murmur2_64(0, key); }
	bool is_empty(size_t key) const { return (key == empty_key); }
	bool is_equal(size_t key1, size_t key2) const { return (key1 == key2); }
};

struct CellInfo {
	uint32_t offset;
	uint32_t idx_count;
	uint32_t vtx_count;
	float max_z;
	float min_z;
	float shift_x;
	float shift_y;
	float shift_z;
};

int init_mgrid(const char *fin, const char *fout)
{
	Mesh mesh;
	MBuf data;
	load_ply(mesh, data, fin);

	constexpr uint32_t DIV = 4;
	TArray<CellInfo> cells(DIV * DIV);
	/* Set cells shift_x and shift_y; */
	for (size_t i = 0; i < DIV * DIV; ++i) {
		cells[i].shift_x = (i & (DIV - 1)) * (1.f / DIV);
		cells[i].shift_y = (i / DIV) * (1.f / DIV);
	}
	/* Init idx_count and vtx_count to 0, and min_z max_z to out of range */
	for (size_t i = 0; i < DIV * DIV; ++i) {
		cells[i].idx_count = 0;
		cells[i].vtx_count = 0;
		cells[i].max_z = -1;
		cells[i].min_z = 9;
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
		int cell_idx = cell_x + cell_y * DIV;
		/* Update max in min altitudes in cell */
		cells[cell_idx].max_z = MAX(cells[cell_idx].max_z, v0.z);
		cells[cell_idx].max_z = MAX(cells[cell_idx].max_z, v1.z);
		cells[cell_idx].max_z = MAX(cells[cell_idx].max_z, v2.z);
		cells[cell_idx].min_z = MIN(cells[cell_idx].min_z, v0.z);
		cells[cell_idx].min_z = MIN(cells[cell_idx].min_z, v1.z);
		cells[cell_idx].min_z = MIN(cells[cell_idx].min_z, v2.z);
		/* Update dico and idx counts */
		tri_idx_to_cell_idx[i] = cell_idx;
		cells[cell_idx].idx_count += 3;
	}
	/* Set shift_z */
	for (size_t i = 1; i < DIV * DIV; ++i) {
		printf("Span : %f\n", cells[i].max_z - cells[i].min_z);
		cells[i].shift_z = floor((1 << 6) * cells[i].min_z) / (1 << 6);
		if (cells[i].shift_z < 0)
			cells[i].shift_z -= 1.f / (1 << 6);
	}
	/* Accumulate idx_counts to offset and reset the latter */
	uint32_t total_idx = 0;
	for (size_t i = 0; i < DIV * DIV; ++i) {
		cells[i].offset = total_idx;
		total_idx += cells[i].idx_count;
		cells[i].idx_count = 0;
	}

	/* Hash table for vertices in cells */
	HashTable<size_t, uint32_t, EncVertexHasher> table(total_idx);
	TArray<uint32_t> cells_idx(total_idx);
	TArray<size_t> cells_vtx(total_idx);
	for (size_t i = 0; i < tri_count; ++i) {
		uint32_t cell_idx = tri_idx_to_cell_idx[i];
		CellInfo &c = cells[cell_idx];
		uint32_t v_idx[3];
		for (int k = 0; k < 3; ++k) {
			uint32_t *pval;
			const Vec3 v = positions[indices[3 * i + k]];
			size_t enc = 0;
			enc += (v.x - c.shift_x) * (1 << 15) * DIV + (1 << 14);
			enc <<= 16;
			enc += (v.y - c.shift_y) * (1 << 15) * DIV + (1 << 14);
			enc <<= 16;
			enc += (v.z - c.shift_z) * (1 << 15);
			enc <<= 16;
			enc += cell_idx;
			pval = table.get_or_set(enc, c.vtx_count);
			size_t *vtx = &cells_vtx[c.offset];
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
		uint32_t *idx = &cells_idx[c.offset];
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
