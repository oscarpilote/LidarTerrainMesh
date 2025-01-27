#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "array.h"
#include "hash.h"
#include "hash_table.h"
#include "math_utils.h"
#include "mesh.h"
#include "mesh_mgrid.h"
#include "mesh_ply.h"
#include "meshoptimizer/src/meshoptimizer.h"
#include "simplifier_mod.h"

#include "swiss_mgrid.h"

struct EncVertex {
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
	static constexpr EncVertex empty_key = {~static_cast<size_t>(0)};
	size_t hash(EncVertex v) const { return murmur2_64(0, v.key); }
	bool is_empty(EncVertex v) const { return (v.key == empty_key.key); }
	bool is_equal(EncVertex v1, EncVertex v2) const
	{
		return (v1.key == v2.key);
	}
};

int init_mgrid(const char *fin, const char *fout, const Pyramid &p)
{
	Mesh mesh;
	MBuf data;
	load_ply(mesh, data, fin);
	const uint32_t DIV = 1 << p.sublevels;
	TArray<CellMeshInfo> info(DIV * DIV);
	TArray<float> shift_x(DIV * DIV);
	TArray<float> shift_y(DIV * DIV);
	TArray<float> max_z(DIV * DIV);
	TArray<float> min_z(DIV * DIV);
	TArray<uint32_t> offset(DIV * DIV);
	/* Initialization */
	for (size_t i = 0; i < DIV * DIV; ++i) {
		info[i].bbox.min.x = 0xffff;
		info[i].bbox.min.y = 0xffff;
		info[i].bbox.min.z = 0xffff;
		info[i].bbox.max.x = 0;
		info[i].bbox.max.y = 0;
		info[i].bbox.max.z = 0;
		info[i].scale_z = 2;
		info[i].vtx_count = 0;
		info[i].idx_count = 0;
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
		// printf("Span : %f\n", max_z[i] - min_z[i]);
		float scal = info[i].scale_z;
		info[i].shift_z = scal * floor(100000 * min_z[i] / scal);
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
	TArray<TVec3<uint16_t>> cells_vtx(total_idx);
	for (size_t i = 0; i < tri_count; ++i) {
		uint32_t cidx = tri_idx_to_cell_idx[i];
		CellMeshInfo &c = info[cidx];
		uint32_t v_idx[3];
		for (int k = 0; k < 3; ++k) {
			uint32_t *pval;
			const Vec3 v = positions[indices[3 * i + k]];
			EncVertex ev;
			ev.ux =
			    (v.x - shift_x[cidx]) * (1 << 15) * DIV + (1 << 14);
			ev.uy =
			    (v.y - shift_y[cidx]) * (1 << 15) * DIV + (1 << 14);
			ev.uz = (100000 * v.z - c.shift_z) / c.scale_z;
			ev.cidx = cidx;
			pval = table.get_or_set(ev, c.vtx_count);
			TVec3<uint16_t> *vtx = &cells_vtx[offset[cidx]];
			if (!pval) {
				v_idx[k] = c.vtx_count;
				vtx[c.vtx_count].x = ev.ux;
				vtx[c.vtx_count].y = ev.uy;
				vtx[c.vtx_count].z = ev.uz;
				c.vtx_count++;
				c.bbox.min.x = MIN(c.bbox.min.x, ev.ux);
				c.bbox.min.y = MIN(c.bbox.min.y, ev.uy);
				c.bbox.min.z = MIN(c.bbox.min.z, ev.uz);
				c.bbox.max.x = MAX(c.bbox.max.x, ev.ux);
				c.bbox.max.y = MAX(c.bbox.max.y, ev.uy);
				c.bbox.max.z = MAX(c.bbox.max.z, ev.uz);
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
	if (!f) {
		printf("Error opening fout\n");
		return -1;
	}
	/* First Pyramid */
	if (!fwrite(&p, sizeof(Pyramid), 1, f)) {
		printf("Error writing pyramid\n");

		fclose(f);
		return -1;
	}
	/* Then CellMeshInfo not yet built */
	CellMeshInfo dummy;
	memset(&dummy, 0, sizeof(CellMeshInfo));
	for (size_t i = 0; i < (DIV * DIV - 1) / 3; ++i) {
		if (!fwrite(&dummy, sizeof(CellMeshInfo), 1, f)) {
			printf("Error writing dummy info\n");
			fclose(f);
			return -1;
		}
	}
	/* Then CellMeshInfo for the just built cells */
	for (size_t i = 0; i < DIV * DIV; ++i) {
		if (!fwrite(&info[i], sizeof(CellMeshInfo), 1, f)) {
			printf("Error writing info\n");
			fclose(f);
			return -1;
		}
	}
	/* Fibally the built cells data */
	for (size_t i = 0; i < DIV * DIV; ++i) {
		TVec3<uint16_t> *vtx = &cells_vtx[offset[i]];
		if (!fwrite(vtx, info[i].vtx_count * sizeof(TVec3<uint16_t>), 1,
			    f)) {
			printf("Error writing vertices\n");
			fclose(f);
			return -1;
		}
		uint32_t *idx = &cells_idx[offset[i]];
		if (!fwrite(idx, info[i].idx_count * sizeof(uint32_t), 1, f)) {
			printf("Error writing indices\n");
			fclose(f);
			return -1;
		}
		/* Write garbage at the moment */
		if (!fwrite(vtx, sizeof(uint32_t) * info[i].vtx_count, 1, f)) {
			printf("Error writing parents %zu\n", i);
			fclose(f);
			return (-1);
		}
	}
	fclose(f);
	printf("Success!\n");
	return (0);
}

/*
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
}*/

/* nx , ny are the coordinates of the block center at level of detail lod */
static int load_blk_cells(MBuf &data, Mesh cells[16], int nx0, int ny0, int lod)
{
	char cur_file[24] = {0};
	char req_file[24] = {0};
	FILE *f = NULL;
	Pyramid p;
	CellMeshInfo info;

	uint32_t blk_idx_count = 0;
	uint32_t blk_vtx_count = 0;

	int cell_count = 0;

	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; j++) {
			int idx_cell = 4 * i + j;
			cells[idx_cell].index_offset = blk_idx_count;
			/* We shall already apply offsets here */
			cells[idx_cell].vertex_offset = 0;
			// cells[idx_cell].vertex_offset = blk_vtx_count;
			int nx = nx0 - 2 * ((i & 1) == 0) + ((j & 1) >> 0);
			int ny = ny0 - 2 * ((i & 2) == 0) + ((j & 2) >> 1);
			int rank = get_mgrid_layout(req_file, p, nx, ny, lod);
			assert(rank >= 0);
			/* Update cur_file if needed */
			if (!cur_file[0] || strcmp(cur_file, req_file) != 0) {
				if (f)
					fclose(f);
				f = fopen(req_file, "rb");
				if (f) {
					strcpy(cur_file, req_file);
				} else {
					cells[idx_cell].index_count = 0;
					cells[idx_cell].vertex_count = 0;
					continue;
				}
			}
			cell_count++;
			/* Load cell */
			/* For vtx we use a small hack to avoid memory alloc */
			if (mg_load_cell(rank, info, NULL, NULL, f) != 0)
				return (1);
			data.reserve_vertices(blk_vtx_count + info.vtx_count);
			uint16_t *vtx =
			    (uint16_t *)(data.positions + blk_vtx_count) +
			    3 * info.vtx_count * sizeof(uint16_t);
			data.reserve_indices(blk_idx_count + info.idx_count);
			uint32_t *idx = data.indices + blk_idx_count;
			if (mg_load_cell(rank, info, vtx, idx, f) != 0)
				return (-1);
			/* Apply iffsets to indices */
			for (uint32_t k = 0; k < info.idx_count; ++k) {
				idx[k] += blk_vtx_count;
			}
			/* Scale/offset vertices */
			int scale_z = info.scale_z;
			int shift_z = info.shift_z;
			Vec3 *blk_vtx = data.positions + blk_vtx_count;
			for (uint32_t k = 0; k < info.vtx_count; ++k) {
				int xi = (int)vtx[3 * k + 0] - (1 << 14) +
					 (nx - nx0) * (1 << 15);
				blk_vtx[k].x = (float)xi / (1 << 16);
				int yi = (int)vtx[3 * k + 1] - (1 << 14) +
					 (ny - ny0) * (1 << 15);
				blk_vtx[k].y = (float)yi / (1 << 16);
				int zi =
				    (int)vtx[3 * k + 2] * scale_z + shift_z;
				blk_vtx[k].z =
				    (float)zi / (100000 * (1 << (lod + 1)));
			}
			cells[idx_cell].index_count = info.idx_count;
			cells[idx_cell].vertex_count = info.vtx_count;
			blk_vtx_count += info.vtx_count;
			blk_idx_count += info.idx_count;
		}
	}
	if (f)
		fclose(f);
	return cell_count;
}

// TODO : This could be a compact mesh vertices utility fct in mesh_utils.cpp
static int join_blk_cells(MBuf &data, uint32_t vtx_count, uint32_t idx_count,
			  uint32_t *remap)
{
	uint32_t unique_vtx_count =
	    meshopt_generateVertexRemap(remap, NULL, vtx_count, data.positions,
					vtx_count, 3 * sizeof(float));
	/* Remap data positions and indices */
	/* Can be safely done in place here. */
	for (uint32_t k = 0; k < vtx_count; ++k) {
		assert(remap[k] <= k);
		data.positions[remap[k]] = data.positions[k];
	}
	for (uint32_t k = 0; k < idx_count; ++k) {
		data.indices[k] = remap[data.indices[k]];
	}
	return unique_vtx_count;
}

static int simplify_block(MBuf &data, uint32_t idx_count, uint32_t vtx_count,
			  float target_error, uint32_t *remap)
{
	TArray<uint32_t> unused(idx_count);
	meshopt_simplify_mod(unused.data, remap, data.indices, idx_count,
			     (const float *)data.positions, vtx_count,
			     3 * sizeof(float), idx_count / 4, target_error,
			     NULL);
	return (0);
}

int build_block(MBuf &data, int nx0, int ny0, int lod)
{
	Mesh cells[16];
	int cell_count = load_blk_cells(data, cells, nx0, ny0, lod);
	uint32_t orig_vtx_count =
	    cells[15].vertex_offset + cells[15].vertex_count;
	uint32_t orig_idx_count =
	    cells[15].index_offset + cells[15].index_count;
	TArray<uint32_t> blk_remap(orig_vtx_count);
	uint32_t blk_vtx_count = join_blk_cells(data, orig_vtx_count,
						orig_idx_count, blk_remap.data);
	float target_error = 0.01;
	TArray<uint32_t> simp_remap(blk_vtx_count);
	simplify_block(data, orig_idx_count, blk_vtx_count, target_error,
		       simp_remap.data);

	return (0);
}

