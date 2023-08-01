#include <assert.h>
#include <stddef.h>
#include <stdint.h>

#include "math_utils.h"
#include "mesh_mgrid.h"

/*
 * The structure of .mg files is :
 *   - MPyr
 *   - MCellInfo[1 + 4 + 4^2 + ... + 4^levels]
 *   - Cells data, in arbitrary order. A cell data contains array of positions,
 *     then array of parents, then array of indices.
 */

struct MCellInfo {
	/* Byte offset from the start of the file */
	size_t offset;
	/* Number of vertices */
	uint32_t vtx_count;
	/* Number of indices */
	uint32_t idx_count;
};

struct EMCellInfo {
	/* Byte offset from the start of the file */
	size_t offset;
	/* Size in bytes of the (encoded) vertex data */
	uint32_t vtx_bufsize;
	/* Number of vertices */
	uint32_t vtx_count;
	/* Size in bytes of the (encoded) index data */
	uint32_t idx_bufsize;
	/* Number of indices */
	uint32_t idx_count;
};

static bool pyramid_contains_cell(MPyr pyr, MCell cell)
{
	bool b = cell.xll >= pyr.base.xll && cell.yll >= pyr.base.yll &&
		 cell.xll + cell.extent <= pyr.base.xll + pyr.base.extent &&
		 cell.yll + cell.extent <= pyr.base.yll + pyr.base.extent &&
		 cell.extent * (1 << pyr.levels) >= pyr.base.extent;
	return b;
}

static int get_pos_in_pyramid(MPyr pyr, MCell cell)
{
	assert(pyramid_contains_cell(pyr, cell));
	int division = pyr.base.extent / cell.extent;
	int rank = (POW2(division) - 1) / 3;
	rank += (cell.xll - pyr.base.xll) / cell.extent;
	rank += (cell.yll - pyr.base.yll) / cell.extent * division;
	return rank;
}

int mg_load_cell(Mesh &m, MBuf &data, MPyr pyr, MCell cell, FILE *f)
{
	int rank = get_pos_in_pyramid(pyr, cell);
	size_t offset;
	offset = sizeof(MPyr) + rank * sizeof(MCellInfo);
	fseek(f, offset, SEEK_SET);
	MCellInfo cinfo;
	if (fread(&cinfo, sizeof(MCellInfo), 1, f) != 1) {
		return -1;
	}
	m.vertex_count = cinfo.vtx_count;
	m.index_count = cinfo.idx_count;
	if (!m.vertex_count || !m.index_count) {
		/* Empty cell is OK*/
		return 0;
	}
	offset = cinfo.offset;

	fseek(f, offset, SEEK_SET);
	void *dest;

	data.reserve_vertices(m.vertex_offset + m.vertex_count);
	dest = data.positions + m.vertex_offset;
	if (fread(dest, sizeof(*data.positions), m.vertex_count, f) !=
	    m.vertex_count) {
		return -1;
	}
	dest = data.remap + m.vertex_offset;
	if (fread(dest, sizeof(*data.remap), m.vertex_count, f) !=
	    m.vertex_count) {
		return -1;
	}

	data.reserve_indices(m.index_offset + m.index_count);
	dest = data.indices + m.index_offset;
	if (fread(dest, sizeof(*data.indices), m.index_count, f) !=
	    m.index_count) {
		return -1;
	}
	return 0;
}

int mg_append_cell(const Mesh &m, const MBuf &data, MPyr pyr, MCell cell,
		   FILE *f)
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

