#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include "mesh.h"
#include "multigrid.h"
//#include "mesh_mgrid.h"

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

struct Cell2DInfo {
	uint32_t vtx_count;
	uint32_t idx_count;
	float shift_z;
	float scale_z;
};

int mg_load_cell(Cell2DInfo &info, void *vtx, void *idx, Cell2D cell,
		 Pyramid2D pyr, Cell2D cell, FILE *f)
{
	int rank = get_cell_index(cell, pyr);
	size_t offset = 0;
	offset += sizeof(Pyramid2D);
	offset += rank * (sizeof(size_t) + sizeof(Cell2DInfo));
	fseek(f, offset, SEEK_SET);
	/* Read data offset */
	if (fread(&offset, sizeof(size_t), 1, f) != 1) {
		return -1;
	}
	/* Read cell info */
	if (fread(&info, sizeof(Cell2DInfo), 1, f) != 1) {
		return -1;
	}
	if (!info.vtx_count || !info.idx_count) {
		/* Empty cell is OK*/
		return 0;
	}
	/* Jump to cell data */
	fseek(f, offset, SEEK_SET);
	/* Read cell data */
	if (!fread(vtx, info.vtx_count * sizeof(uint16_t), 1, f)) {
		return -1;
	}
	if (!fread(idx, info.idx_count * sizeof(uint32_t), 1, f)) {
		return -1;
	}
	return 0;
}

int mg_append_cell(const Cell2D &cellMesh &m, const MBuf &data, MPyr pyr,
		   MCell cell, FILE *f)
{

	fseek(f, 0, SEEK_END);
	size_t offset = ftell(f);
	void *orig;

	orig = data.positions + m.vertex_offset;
	if (fwrite(orig, sizeof(*data.positions), m.vertex_count, f) !=
	    m.vertex_count) {
		return -1;
	}
	orig = data.remap + m.vertex_offset;
	if (fwrite(orig, sizeof(*data.remap), m.vertex_count, f) !=
	    m.vertex_count) {
		return -1;
	}

	orig = data.indices + m.index_offset;
	if (fwrite(orig, sizeof(*data.indices), m.index_count, f) !=
	    m.index_count) {
		return -1;
	}

	MCellInfo cinfo;
	cinfo.vtx_count = m.vertex_count;
	cinfo.idx_count = m.index_count;
	cinfo.offset = offset;

	int rank = get_pos_in_pyramid(pyr, cell);
	offset = sizeof(MPyr) + rank * sizeof(MCellInfo);
	fseek(f, 0, SEEK_SET);
	if (fwrite(&cinfo, sizeof(MCellInfo), 1, f) != 1) {
		return -1;
	}
	return 0;
}

