#include <memory>
#include <stdint.h>
#include <stdio.h>

#include "aabb.h"

struct CopcReader;

struct CopcReader *copc_init(const char *filename);

void copc_fini(CopcReader *copc);

uint32_t copc_set_target_bbox(CopcReader *copc, const TAabb<double> &box);

int copc_cell_point_count(CopcReader *copc, uint32_t cell_idx);

int copc_read_cell(CopcReader *copc, uint32_t cell_idx, char *dst);

uint32_t copc_bound_inside(CopcReader *copc, const TAabb<double> &box);

uint32_t copc_load_inside(CopcReader *copc, const TAabb<double> &box,
			  char *buf);

