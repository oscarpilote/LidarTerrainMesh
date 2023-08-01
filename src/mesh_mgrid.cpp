#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include "mesh_mgrid.h"

/*
 * The structure of .mg files is :
 *   - Pyramid2D
 *   - struct { size_t offset; Cell2DInfo info}[1 + 4 + 4^2 + ... + 4^levels]
 *   - Cells data, in arbitrary order. A cell data contains array of uint16_t
 *     positions, then an array of uint32_t parents, then array of uint32_t
 *     indices.
 */

/*
 * The structure of .emg files is :
 *   - Pyramid2D
 *   - struct {Offset, vtx_bufsize, idx_bufside}[1 + 4 + 4^2 + ... + 4^levels]
 *   - Cell2DInfo[1 + 4 + 4^2 + ... + 4^levels]
 *   - Cells data, in arbitrary order. A cell data contains an encoding of
 *     vertex data followed by an encoding of index data.
 */

int mg_load_cell(unsigned rank, CellMeshInfo &info, void *vtx, void *idx,
		 FILE *f)
{
	size_t meta_offset =
	    sizeof(Pyramid) + rank * (sizeof(size_t) + sizeof(CellMeshInfo));
	fseek(f, meta_offset, SEEK_SET);
	/* Read data offset */
	size_t data_offset;
	if (!fread(&data_offset, sizeof(size_t), 1, f)) {
		return -1;
	}
	/* Read cell info */
	if (!fread(&info, sizeof(CellMeshInfo), 1, f)) {
		return -1;
	}
	if (!info.vtx_count || !info.idx_count) {
		/* Empty cell is OK*/
		return 0;
	}
	/* If vtx or idx is NULL exit early (info only) */
	if (vtx == NULL || idx == NULL) {
		return 0;
	}
	/* Jump to cell data */
	fseek(f, data_offset, SEEK_SET);
	/* Read cell data */
	if (!fread(vtx, 3 * info.vtx_count * sizeof(uint16_t), 1, f)) {
		return -1;
	}
	if (!fread(idx, info.idx_count * sizeof(uint32_t), 1, f)) {
		return -1;
	}
	return 0;
}

int mg_append_cell(unsigned rank, const CellMeshInfo &info, const uint16_t *vtx,
		   const uint32_t *idx, const uint32_t *remap, FILE *f)

{
	/* Get file size */
	fseek(f, 0, SEEK_END);
	size_t data_offset = ftell(f);

	/* Append data at end of file */
	if (!fwrite(vtx, 3 * info.vtx_count * sizeof(uint16_t), 1, f)) {
		return -1;
	}
	if (!fwrite(idx, info.idx_count * sizeof(uint32_t), 1, f)) {
		return -1;
	}
	if (remap) {
		if (!fwrite(remap, info.vtx_count * sizeof(uint32_t), 1, f)) {
			return -1;
		}
	} else {
		uint32_t dummy = (-1);
		uint32_t count = info.vtx_count;
		if (fwrite(&dummy, sizeof(uint32_t), count, f) != count) {
			return -1;
		}
	}
	/* Rewind to write metadata */
	size_t meta_offset =
	    sizeof(Pyramid) + rank * (sizeof(size_t) + sizeof(CellMeshInfo));
	fseek(f, meta_offset, SEEK_SET);
	/* Write data offset */
	if (!fwrite(&data_offset, sizeof(size_t), 1, f)) {
		return -1;
	}
	if (!fwrite(&info, sizeof(CellMeshInfo), 1, f)) {
		return -1;
	}
	return 0;
}

int mg_update_cell(unsigned rank, const uint32_t *remap, FILE *f)

{
	/* Read cell metadata */
	size_t meta_offset =
	    sizeof(Pyramid) + rank * (sizeof(size_t) + sizeof(CellMeshInfo));
	fseek(f, meta_offset, SEEK_SET);
	/* Read data offset */
	size_t data_offset;
	if (!fread(&data_offset, sizeof(size_t), 1, f)) {
		return -1;
	}
	/* Read cell info */
	CellMeshInfo info;
	if (!fread(&info, sizeof(CellMeshInfo), 1, f)) {
		return -1;
	}
	if (!info.vtx_count || !info.idx_count) {
		/* Empty cell is OK*/
		return 0;
	}
	/* Jump to cell data */
	size_t remap_offset = data_offset +
			      3 * info.vtx_count * sizeof(uint16_t) +
			      info.idx_count * sizeof(uint32_t);
	fseek(f, remap_offset, SEEK_SET);
	if (!fwrite(remap, info.vtx_count * sizeof(uint32_t), 1, f)) {
		return -1;
	}
	return 0;
}

