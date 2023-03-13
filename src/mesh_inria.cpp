#include <stdio.h>
#include <string.h>

#include "mmg/mmgs/libmmgs.h"

#include "fast_parse.h"

#include "mesh_inria.h"

#define LINE 64

static inline bool line_should_be_skipped(const char *line)
{
	return (line[0] == '\n' || line[0] == '#' || line[0] == ' ');
}

static int read_positions(char *line, FILE *f, Mesh &mesh, MBuf &data)
{
	fgets(line, LINE, f);
	if (sscanf(line, "%d", &mesh.vertex_count) != 1) {
		printf("Could not read vertex count.\n");
		return (-1);
	}
	data.update_vtx_attr(data.vtx_attr | VtxAttr::POS);
	data.reserve_vertices(mesh.vertex_count);
	const char *ptr;
	for (size_t i = 0; i < mesh.vertex_count; ++i) {
		fgets(line, LINE, f);
		ptr = line;
		ptr = parse_float(ptr, &data.positions[i].x);
		ptr = parse_float(ptr, &data.positions[i].y);
		ptr = parse_float(ptr, &data.positions[i].z);
	}
	return (0);
}

static int read_indices(char *line, FILE *f, Mesh &mesh, MBuf &data)
{
	fgets(line, LINE, f);
	unsigned tri_count;
	if (sscanf(line, "%d", &tri_count) != 1) {
		printf("Could not read triangle count.\n");
		return (-1);
	}
	mesh.index_count = 3 * tri_count;
	data.reserve_indices(mesh.index_count);
	const char *ptr;
	for (size_t i = 0; i < tri_count; ++i) {
		fgets(line, LINE, f);
		ptr = line;
		ptr = parse_int(ptr, (int *)&data.indices[3 * i + 0]);
		ptr = parse_int(ptr, (int *)&data.indices[3 * i + 1]);
		ptr = parse_int(ptr, (int *)&data.indices[3 * i + 2]);
		--data.indices[3 * i + 0];
		--data.indices[3 * i + 1];
		--data.indices[3 * i + 2];
	}
	return (0);
}

static int skip_section(char *line, FILE *f)
{
	fgets(line, LINE, f);
	unsigned dummy_count;
	if (sscanf(line, "%d", &dummy_count) != 1) {
		printf("Could not read dummy count : %s\n", line);
		return (-1);
	}
	for (size_t i = 0; i < dummy_count; ++i) {
		fgets(line, LINE, f);
	}
	return (0);
}

int read_inria(Mesh &mesh, MBuf &data, const char *fname)
{

	// TODO : rethink about this clear
	mesh.clear();
	data.clear();

	FILE *f = fopen(fname, "r");
	if (f == NULL)
		return (-1);

	char line[LINE];

	float ver;
	int dim;

	int section = 0;

	while (fgets(line, LINE, f)) {
		if (line_should_be_skipped(line))
			continue;
		switch (section) {
		case 0:
			if (sscanf(line, "MeshVersionFormatted %f", &ver) !=
			    1) {
				printf("Not a valid .mesh file : %s\n", fname);
				fclose(f);
				return (-1);
			}
			section = 1;
			break;
		case 1:
			if (sscanf(line, "Dimension %d", &dim) != 1) {
				printf("Not a valid .mesh file "
				       ": %s\n",
				       fname);
				fclose(f);
				return (-1);
			}
			if (dim != 3) {
				printf("Only dimension 3 is "
				       "supported by this "
				       "reader.\n");
				fclose(f);
				return (-1);
			}
			section = 2;
			break;
		case 2:
			if (strncmp(line, "Vertices", 8) == 0) {
				read_positions(line, f, mesh, data);
			} else if (strncmp(line, "Triangles", 9) == 0) {
				read_indices(line, f, mesh, data);
			} else if (strncmp(line, "End", 3) == 0) {
				fclose(f);
				return (0);
			} else {
				skip_section(line, f);
			}
		}
	}
	fclose(f);
	return (0);
}
#undef LINE

int write_inria(const char *fname, const Mesh &mesh, const MBuf &data)
{
	FILE *f = fopen(fname, "w");
	if (f == NULL)
		return (-1);
	fprintf(f, "MeshVersionFormatted 1\n");
	fprintf(f, "Dimension 3\n\n");
	fprintf(f, "Vertices\n");
	fprintf(f, "%d\n", mesh.vertex_count);
	for (size_t i = 0; i < mesh.vertex_count; ++i) {
		fprintf(f, "%.12g %.12g %.12g 0\n", data.positions[i].x,
			data.positions[i].y, data.positions[i].z);
	}
	if (data.vtx_attr & VtxAttr::NML) {
		fprintf(f, "\n");
		fprintf(f, "Normals\n");
		fprintf(f, "%d\n", mesh.vertex_count);
		for (size_t i = 0; i < mesh.vertex_count; ++i) {
			fprintf(f, "%.12g %.12g %.12g\n", data.normals[i].x,
				data.normals[i].y, data.normals[i].z);
		}
	}
	fprintf(f, "\n");
	fprintf(f, "Triangles\n");
	unsigned tri_count = mesh.index_count / 3;
	fprintf(f, "%d\n", tri_count);
	for (size_t i = 0; i < tri_count; ++i) {
		fprintf(f, "%d %d %d %d\n", data.indices[3 * i + 0] + 1,
			data.indices[3 * i + 1] + 1,
			data.indices[3 * i + 2] + 1, 0);
	}
	fclose(f);
	return 0;
}

int load_mesh_to_mmg(const Mesh &mesh, const MBuf &data, MMG5_pMesh mm,
		     MMG5_pSol ss)
{
	MMGS_Set_meshSize(mm, mesh.vertex_count, mesh.index_count / 3, 0);
	const Vec3 *pos = data.positions + mesh.vertex_offset;
	for (size_t i = 0; i < mesh.vertex_count; ++i) {
		double c0 = pos[i].x;
		double c1 = pos[i].y;
		double c2 = pos[i].z;
		MMGS_Set_vertex(mm, c0, c1, c2, 0, i + 1);
	}
	const uint32_t *idx = data.indices + mesh.index_offset;
	for (size_t i = 0; i < mesh.index_count / 3; ++i) {
		uint32_t i0 = idx[3 * i + 0];
		uint32_t i1 = idx[3 * i + 1];
		uint32_t i2 = idx[3 * i + 2];
		MMGS_Set_triangle(mm, i0 + 1, i1 + 1, i2 + 1, 0, i + 1);
	}

	if (MMGS_Chk_meshData(mm, ss) != 1)
		return (-1);

	return 0;
}

int unload_mesh_from_mmg(Mesh &mesh, MBuf &data, const MMG5_pMesh mm,
			 const MMG5_pSol ss)
{
	int np, nt, na;
	MMGS_Get_meshSize(mm, &np, &nt, &na);

	data.reserve_vertices(mesh.vertex_offset + np);
	Vec3 *pos = data.positions + mesh.vertex_offset;
	for (int i = 0; i < np; ++i) {
		double c0, c1, c2;
		int dummy;
		MMGS_GetByIdx_vertex(mm, &c0, &c1, &c2, &dummy, &dummy, &dummy,
				     i + 1);
		pos[i].x = c0;
		pos[i].y = c1;
		pos[i].z = c2;
	}
	mesh.vertex_count = np;

	data.reserve_indices(mesh.index_offset + 3 * nt);
	uint32_t *idx = data.indices + mesh.index_offset;
	for (int i = 0; i < nt; ++i) {
		int i0, i1, i2;
		int dummy;
		MMGS_Get_triangle(mm, &i0, &i1, &i2, &dummy, &dummy);
		idx[3 * i + 0] = i0 - 1;
		idx[3 * i + 1] = i1 - 1;
		idx[3 * i + 2] = i2 - 1;
	}
	mesh.index_count = 3 * nt;

	return 0;
}

int mmg_remesh(Mesh &mesh, MBuf &data, float hausd, float hgrad, bool ridges,
	       int verbose)
{
	MMG5_pMesh mm = NULL;
	MMG5_pSol ss = NULL;
	MMGS_Init_mesh(MMG5_ARG_start, MMG5_ARG_ppMesh, &mm, MMG5_ARG_ppMet,
		       &ss, MMG5_ARG_end);

	MMGS_Init_parameters(mm);

	MMGS_Set_dparameter(mm, ss, MMGS_DPARAM_hausd, hausd);
	MMGS_Set_dparameter(mm, ss, MMGS_DPARAM_hgrad, hgrad);
	MMGS_Set_iparameter(mm, ss, MMGS_IPARAM_angle, ridges);
	MMGS_Set_iparameter(mm, ss, MMGS_IPARAM_verbose, verbose);

	load_mesh_to_mmg(mesh, data, mm, ss);

	MMGS_mmgslib(mm, ss);
	int np, nt, na;
	MMGS_Get_meshSize(mm, &np, &nt, &na);

	unload_mesh_from_mmg(mesh, data, mm, ss);

	MMGS_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh, &mm, MMG5_ARG_ppMet, &ss,
		      MMG5_ARG_end);

	return (0);
}
