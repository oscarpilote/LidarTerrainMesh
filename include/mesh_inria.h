#pragma once

#include "mesh.h"



int read_inria(Mesh &mesh, MBuf &data, const char *fname);
int write_inria(const char *fname, const Mesh &mesh, const MBuf &data);

int mmg_remesh(Mesh &m, MBuf &data, float hausd, float hgrad = 5.f, bool ridges = false, int verbose = -1);

