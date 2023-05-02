#pragma once
#include "multigrid.h"

int init_mgrid(const char *fin, const char *fout, const Pyramid &p);

int get_mgrid_layout(char filename[24], Pyramid &p, int nx, int ny, int zl);
