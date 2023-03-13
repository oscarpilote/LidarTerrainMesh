#pragma once

#include "mmg/mmgs/libmmgs.h"

#include "mesh.h"

int read_inria(Mesh &mesh, MBuf &data, const char *fname);
int write_inria(const char *fname, const Mesh &mesh, const MBuf &data);

int load_mesh_to_mmg(const Mesh &m, const MBuf &data, MMG5_pMesh mm,
		     MMG5_pSol ss);
int unload_mesh_from_mmg(Mesh &m, MBuf &data, const MMG5_pMesh mm,
			 const MMG5_pSol ss);
int mmg_remesh(Mesh &m, MBuf &data, float hausd, float hgrad = 5.f,
	       bool ridges = false, int verbose = -1);

