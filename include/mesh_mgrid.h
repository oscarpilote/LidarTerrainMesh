#pragma once

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include "mesh.h"

struct MCell {
	float xll;
	float yll;
	float extent;
};

struct MPyr {
	MCell base;
	int levels;
};

int mg_load_cell(Mesh &m, MBuf &data, MPyr pyr, MCell cell, FILE *f);

int mg_append_cell(const Mesh &m, const MBuf &data, MPyr pyr, MCell cell,
		   FILE *f);

