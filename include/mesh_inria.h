#pragma once

#include "mesh.h"

int read_inria(Mesh &mesh, MBuf &data, const char *fname);
int write_inria(const char *fname, const Mesh &mesh, const MBuf &data);

