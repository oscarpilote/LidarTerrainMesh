#pragma once

#include "mesh.h"

int load_ply(Mesh &mesh, MBuf &data, const char *fname,
	     int filter = VtxAttr::ALL);
int write_ply(const char *fname, const Mesh &mesh, const MBuf &data,
	      int filter = VtxAttr::ALL);

