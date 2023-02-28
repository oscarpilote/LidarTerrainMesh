#include <stdio.h>
#include <string.h>

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
		printf("Could not read dummy count.\n");
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
		fprintf(f, "%lf %lf %lf 0\n", data.positions[i].x,
			data.positions[i].y, data.positions[i].z);
	}
	if (data.vtx_attr & VtxAttr::NML) {
		fprintf(f, "\n");
		fprintf(f, "Normals\n");
		fprintf(f, "%d\n", mesh.vertex_count);
		for (size_t i = 0; i < mesh.vertex_count; ++i) {
			fprintf(f, "%lf %lf %lf\n", data.normals[i].x,
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
