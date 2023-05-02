#pragma once

#include <stdint.h>

struct CellCoord {
	int nx;
	int ny;
	int lod;
};

struct Pyramid {
	CellCoord top;
	int sublevels;
};

inline struct CellCoord cell_child(struct CellCoord parent, int k)
{
	struct CellCoord child;
	child.nx = 2 * parent.nx + (k & 1);
	child.ny = 2 * parent.ny + (k & 2 ? 1 : 0);
	child.lod = parent.lod - 1;
	return child;
}

inline struct CellCoord cell_parent(struct CellCoord child)
{
	struct CellCoord parent;
	parent.nx = (child.nx - (child.nx < 0)) / 2;
	parent.ny = (child.ny - (child.ny < 0)) / 2;
	parent.lod = child.lod + 1;
	return parent;
}

inline int get_cell_index(const CellCoord &c, const Pyramid &p)
{
	if (c.lod > p.top.lod)
		return (-1);
	if (c.lod < p.top.lod - p.sublevels) {
		return (-1);
	}
	int N = 1 << (p.top.lod - c.lod);
	int di = c.nx - (p.top.nx * N);
	int dj = c.ny - (p.top.ny * N);
	if (di < 0 || di >= N || dj < 0 || dj >= N) {
		return (-1);
	}
	int skip = (N * N - 1) / 3;
	return (skip + dj * N + di);
}
