#include <stddef.h>
#include <stdint.h>

#include "mesh_mgrid.h"

struct CellInfo {
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

