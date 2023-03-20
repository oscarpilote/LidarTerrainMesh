#include "kd_tree.h"
#include <stdio.h>

#include "mesh.h"
#include "mesh_adjacency.h"
#include "mesh_ply.h"
#include "vec3.h"

int main(int argc, char **argv)
{

	if (argc < 5) {
		return -1;
	}

	Mesh mesh[2];
	MBuf data[2];
	EAdj adj[2];

	TArray<bool> is_bdy[2];
	TArray<uint32_t> candidates[2];
	uint32_t counts[2];
	for (int i = 0; i < 2; ++i) {
		load_ply(mesh[i], data[i], argv[i + 1]);
		fill_edge_adjacency(mesh[i], data[i], adj[i]);
		counts[i] = find_boundary(adj[i], is_bdy[i]);
		const Vec3 *pos = data[i].positions + mesh[i].vertex_offset;
		candidates[i].reserve(counts[i] / 4 + 32);
		for (size_t j = 0; j < adj[i].vertex_count; ++j) {
			if (is_bdy[i][j] && (pos[j].x == (1 - i))) {
				candidates[i].push_back(j);
			}
		}
		counts[i] = candidates[i].size;
	}

	TArray<Vec3> kd_pos[2];
	for (int i = 0; i < 2; ++i) {
		const Vec3 *pos = data[i].positions + mesh[i].vertex_offset;
		kd_pos[i].resize(counts[i]);
		for (size_t j = 0; j < counts[i]; ++j) {
			kd_pos[i][j] = pos[candidates[i][j]];
		}
	}
	KdCoords<3> kdcoords0{(const float *)&kd_pos[0][0], counts[0]};
	KdCoords<3> kdcoords1{(const float *)&kd_pos[1][0], counts[1]};
	KdTree<3> tree0(3, kdcoords0, 10);
	KdTree<3> tree1(3, kdcoords1, 10);

	/* Stitch mesh1 vertices to mesh0 one */
	uint32_t knn_idx;
	float sqr_dist;
	TArray<uint32_t> remap1(counts[1]);
	TArray<bool> stitched(counts[0], false);
	for (size_t j = 0; j < counts[1]; ++j) {
		const float *query = &kd_pos[1][j].x;
		tree0.knnSearch(query, 1, &knn_idx, &sqr_dist);
		remap1[j] = knn_idx;
		stitched[knn_idx] = true;
	}
	/* Stitch mesh0 vertices that are not alread done */
	TArray<uint32_t> remap0(counts[0]);
	for (size_t j = 0; j < counts[0]; ++j) {
		if (stitched[j]) {
			remap0[j] = j;
		} else {
			const float *query = &kd_pos[0][j].x;
			tree1.knnSearch(query, 1, &knn_idx, &sqr_dist);
			remap0[j] = remap1[knn_idx];
		}
	}

	/* Modify initial  positions in to stitched ones */
	Vec3 *pos0 = data[0].positions + mesh[0].vertex_offset;
	Vec3 *pos1 = data[1].positions + mesh[1].vertex_offset;
	for (size_t j = 0; j < counts[0]; ++j) {
		if (remap0[j] != j) {
			pos0[candidates[0][j]] = pos0[candidates[0][remap0[j]]];
		}
	}
	for (size_t j = 0; j < counts[1]; ++j) {
		pos1[candidates[1][j]].y = pos0[candidates[0][remap1[j]]].y;
		pos1[candidates[1][j]].z = pos0[candidates[0][remap1[j]]].z;
	}

	/* TODO : optimize to remove duplicated vertices */
	write_ply(argv[3], mesh[0], data[0]);
	write_ply(argv[4], mesh[1], data[1]);

	printf("Done !\n");
	return (0);
}
