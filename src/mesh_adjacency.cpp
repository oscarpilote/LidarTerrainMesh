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

void fill_edge_adjacency(const Mesh &mesh, const MBuf &data, EAdj &adj)
{
	adj.vertex_count = mesh.vertex_count;
	adj.index_count = mesh.index_count;
	adj.edges.resize(mesh.vertex_count + mesh.index_count);
	for (size_t i = 0; i < mesh.vertex_count + mesh.index_count; ++i) {
		adj.edges[i].dest = adj.edges[i].idx = ~0;
	}

	const uint32_t *indices = data.indices + mesh.index_offset;
	const uint32_t off = adj.vertex_count;
	TArray<Edge> &edges = adj.edges;

	size_t tri_count = mesh.index_count / 3;
	for (size_t i = 0; i < tri_count; ++i) {

		uint32_t i0 = indices[3 * i + 0];
		uint32_t i1 = indices[3 * i + 1];
		uint32_t i2 = indices[3 * i + 2];

		uint32_t dest, idx;

		dest = edges[i0].dest;
		idx = edges[i0].idx;
		edges[i0].dest = i1;
		edges[i0].idx = 3 * i + 0;
		edges[off + 3 * i + 0].dest = dest;
		edges[off + 3 * i + 0].idx = idx;

		dest = edges[i1].dest;
		idx = edges[i1].idx;
		edges[i1].dest = i2;
		edges[i1].idx = 3 * i + 1;
		edges[off + 3 * i + 1].dest = dest;
		edges[off + 3 * i + 1].idx = idx;

		dest = edges[i2].dest;
		idx = edges[i2].idx;
		edges[i2].dest = i0;
		edges[i2].idx = 3 * i + 2;
		edges[off + 3 * i + 2].dest = dest;
		edges[off + 3 * i + 2].idx = idx;
	}
}

uint32_t find_edge(const EAdj &adj, uint32_t i0, uint32_t i1)
{
	assert(i0 < adj.vertex_count);

	const TArray<Edge> &edges = adj.edges;
	const uint32_t off = adj.vertex_count;

	uint32_t dest = edges[i0].dest;
	uint32_t idx = edges[i0].idx;

	while (dest != ~0) {
		if (dest == i1) {
			return idx;
		}
		dest = edges[off + idx].dest;
		idx = edges[off + idx].idx;
	}

	return (~0);
}

bool has_edge(const EAdj &adj, uint32_t i0, uint32_t i1)
{
	return (find_edge(adj, i0, i1) != ~0);
}

uint32_t next_edge(const EAdj &adj, uint32_t edge_idx, uint32_t i1)
{
	assert(edge_idx < adj.index_count);

	const TArray<Edge> &edges = adj.edges;
	const uint32_t off = adj.vertex_count;

	uint32_t dest = edges[off + edge_idx].dest;
	uint32_t idx = edges[off + edge_idx].idx;

	while (dest != ~0) {
		if (dest == i1) {
			return idx;
		}
		dest = edges[off + idx].dest;
		idx = edges[off + idx].idx;
	}

	return (~0);
}

uint32_t

    uint32_t
    classify_edges(const Mesh &mesh, const MBuf &data, const EAdj &adj)
{
}

void fill_tri_adjacency(const Mesh &mesh, const MBuf &data, const EAdj eadj,
			TArray<uint32_t> triadj)
{
	triadj.resize(mesh.index_count);
	for (size_t i = 0; i < mesh.index_count; ++i) {
		triadj[i] = ~0;
	}

	const uint32_t *indices = data.indices + mesh.index_offset;
	size_t tri_count = mesh.index_count / 3;
	for (size_t i = 0; i < tri_count; ++i) {
		uint32_t i0 = indices[3 * i + 0];
		uint32_t i1 = indices[3 * i + 1];
		uint32_t i2 = indices[3 * i + 2];
		triadj[3 * i + 0] = find_edge(eadj, i1, i0);
		triadj[3 * i + 1] = find_edge(eadj, i2, i1);
		triadj[3 * i + 2] = find_edge(eadj, i0, i2);
	}
}

uint32_t find_connected_components(const TArray<uint32_t> &triadj,
				   TArray<uint32_t> &cc)
{
	size_t tri_count = triadj.size / 3;

	cc.resize(tri_count);
	for (size_t i = 0; i < tri_count; ++i) {
		cc[i] = ~0;
	}

	uint32_t classified = 0;
	/* FIFO queue */
	TArray<bool> visited(tri_count, false);
	TArray<uint32_t> to_visit(tri_count);

	uint32_t next_tri = 0;
	uint32_t next_cc = 0;
	while (classified != tri_count) {
		uint32_t front = 0;
		uint32_t back = 0;
		while
			uint32_t t1 = triadj[3 * next_tri + 0] / 3;
		uint32_t t2 = triadj[3 * next_tri + 1] / 3;
		uint32_t t3 = triadj[3 * next_tri + 2] / 3;
		if (t1 != ~0)
	}

