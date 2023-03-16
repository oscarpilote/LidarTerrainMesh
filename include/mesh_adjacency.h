#pragma once

#include "array.h"
#include "mesh.h"

struct Edge {
	uint32_t dest;
	uint32_t idx;
};

struct EAdj {
	uint32_t vertex_count;
	uint32_t index_count;
	TArray<Edge> edges;
};

void fill_edge_adjacency(const Mesh &mesh, const MBuf &data, EAdj &adj);

uint32_t find_edge(const EAdj &adj, uint32_t i0, uint32_t i1);

bool has_edge(const EAdj &adj, uint32_t i0, uint32_t i1);

uint32_t next_edge(const EAdj &adj, uint32_t edge_idx, uint32_t i1);

void fill_tri_adjacency(const Mesh &mesh, const MBuf &data, const EAdj &eadj,
			TArray<uint32_t> &triadj);

uint32_t find_connected_components(const TArray<uint32_t> &triadj,
				   TArray<uint32_t> &cc);

