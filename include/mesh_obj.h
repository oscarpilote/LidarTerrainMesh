#pragma once

#include "mesh.h"

int load_obj(Mesh &mesh, MBuf &data, const char *fname);
int write_obj(const Mesh &m, const MBuf &data, const char *fname);
