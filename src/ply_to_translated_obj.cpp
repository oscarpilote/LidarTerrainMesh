#include <stdio.h>

#include "mesh.h"
#include "mesh_obj.h"
#include "mesh_ply.h"

int main(int argc, char **argv)
{
	if (argc < 5) {
		printf("Syntax : %s src_ply dst_obj dx dy\n", argv[0]);
		return -1;
	}
	char *src = argv[1];
	char *dst = argv[2];
	float dx = atof(argv[3]);
	float dy = atof(argv[4]);
	Mesh mesh;
	MBuf data;
	load_ply(mesh, data, src);
	for (size_t i = 0; i < mesh.vertex_count; i++) {
		data.positions[i].x += dx;
		data.positions[i].y += dy;
	}
	write_obj(mesh, data, dst);
	return 0;
}
