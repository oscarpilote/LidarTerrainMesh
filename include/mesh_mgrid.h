#pragma once

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include "multigrid.h"
#include "aabb.h"

struct CellMeshInfo {
	/* Bounds in uin16_t encoded x,y,z  */ 
	TAabb<uint16_t> bbox;
	/* Offset and scale in z, in centimeters */
	int shift_z;
	int scale_z;
	/* Vertex and index counts */
	uint32_t vtx_count;
	uint32_t idx_count;
};

int mg_load_cell(unsigned rank, CellMeshInfo &info, void *vtx, void *idx, FILE *f);
int mg_append_cell(unsigned rank, const CellMeshInfo &info, const uint16_t *vtx,
		   const uint32_t *idx, const uint32_t *remap, FILE *f);
int mg_update_cell(unsigned rank, const uint32_t *remap, FILE *f);

