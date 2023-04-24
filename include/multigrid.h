#pragma once

#include <stdint.h>

struct Cell2D {
	int nx;
	int ny;
	int zl;
};

struct Pyramid2D {
	Cell2D base;
	int sublevels;
};

inline struct Cell2D cell_child(struct Cell2D parent, int k)
{
	struct Cell2D child;
	child.nx = 2 * parent.nx + (k & 1);
	child.ny = 2 * parent.ny + (k & 2 ? 1 : 0);
	child.zl = parent.zl + 1;
	return child;
}

inline struct Cell2D cell_parent(struct Cell2D child)
{
	struct Cell2D parent;
	parent.nx = (child.nx - (child.nx < 0)) / 2;
	parent.ny = (child.ny - (child.ny < 0)) / 2;
	parent.zl = child.zl - 1;
	return parent;
}

int get_cell_index(const Cell2D &cell, const Pyramid2D p)
{
	if (cell.zl < p.base.zl)
		return (-1);
	if (cell.zl > p.base.zl + p.sublevels) {
		return (-1);
	}
	int N = 1 << (cell.zl - p.base.zl);
	int di = cell.nx - (p.base.nx * N);
	int dj = cell.ny - (p.base.ny * N);
	if (di < 0 || di >= N || dj < 0 || dj >= N) {
		return (-1);
	}
	int skip = (N * N - 1) / 3;
	return (skip + dj * N + di);
}
