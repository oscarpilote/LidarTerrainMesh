#pragma once

#include <stddef.h>
#include <stdint.h>

#include "mesh.h"
#include "multigrid.h"

int load_pyramid(MPyramid &mp, MBuf &data, const char *fname);
int append int save_pyramid(const MPyramid &mp, const MBuf &data,
			    const char *fname);

