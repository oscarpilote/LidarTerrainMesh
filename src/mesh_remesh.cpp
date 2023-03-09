#include <stdio.h>


#include "mesh.h"
#include "mesh_inria.h"

int mesh_remesh(Mesh& m, MBuf &data, float hausd, float hgrad, bool ridges, int verbose)
{
	return mmg_remesh(m, data, hausd, hgrad, ridges, verbose ? 1 : -1);
}
