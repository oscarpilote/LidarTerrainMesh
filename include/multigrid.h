#pragma once

#include <stdint.h>

struct CellCoord {
	int nx;
	int ny;
	int zl;
};

struct Pyramid {
	CellCoord base;
	int sublevels;
};

inline struct CellCoord cell_child(struct CellCoord parent, int k)
{
	struct CellCoord child;
	child.nx = 2 * parent.nx + (k & 1);
	child.ny = 2 * parent.ny + (k & 2 ? 1 : 0);
	child.zl = parent.zl + 1;
	return child;
}

inline struct CellCoord cell_parent(struct CellCoord child)
{
	struct CellCoord parent;
	parent.nx = (child.nx - (child.nx < 0)) / 2;
	parent.ny = (child.ny - (child.ny < 0)) / 2;
	parent.zl = child.zl - 1;
	return parent;
}

inline int get_cell_index(const CellCoord &c, const Pyramid &p)
{
	if (c.zl < p.base.zl)
		return (-1);
	if (c.zl > p.base.zl + p.sublevels) {
		return (-1);
	}
	int N = 1 << (c.zl - p.base.zl);
	int di = c.nx - (p.base.nx * N);
	int dj = c.ny - (p.base.ny * N);
	if (di < 0 || di >= N || dj < 0 || dj >= N) {
		return (-1);
	}
	int skip = (N * N - 1) / 3;
	return (skip + dj * N + di);
}
