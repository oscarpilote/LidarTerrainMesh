#pragma once

#include <stdint.h>

struct Cell2D {
	int i;
	int j;
	int level;
};

inline struct Cell2D cell_child(struct Cell2D parent, int k)
{
	struct Cell2D child;
	child.i = 2 * parent.i + (k & 1);
	child.j = 2 * parent.j + (k & 2 ? 1 : 0);
	child.level = parent.level - 1;
	return child;
}

inline struct Cell2D cell_parent(struct Cell2D child)
{
	struct Cell2D parent;
	parent.i = (child.i - (child.i < 0)) / 2;
	parent.j = (child.j - (child.j < 0)) / 2;
	parent.level = child.level + 1;
	return parent;
}

int get_subcell_index(struct Cell2D cell, struct Cell2D subcell)
{
	int delta = cell.level - subcell.level;
	if (delta < 0)
		return (-1);
	int di = subcell.i - (cell.i << delta);
	int dj = subcell.j - (cell.j << delta);
	int N = (1 << delta);
	if (di < 0 || di >= N || dj < 0 || dj >= N)
		return (-1);
	int skip = ((1 << (2 * delta)) - 1) / 3;
	return (skip + dj * N + di);
}
