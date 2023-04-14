#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <pthread.h>

#include "meshoptimizer/src/meshoptimizer_mod.h"

#include "camera.h"
#include "hash_table.h"
#include "math_utils.h"
#include "mesh.h"
#include "mesh_grid.h"
#include "mesh_utils.h"
#include "vec3.h"

#include "chrono.h"

static inline void point_to_cell_coord(CellCoord &coord, const Vec3 &p,
				       const Vec3 &base, float inv_step)
{
	Vec3 my_test = (p - base) * inv_step;

	coord.x = floor(my_test.x);
	coord.y = floor(my_test.y);
	coord.z = floor(my_test.z);
}

CellCoord parent_coord(const CellCoord coord)
{
	CellCoord pcoord;

	pcoord.lod = coord.lod + 1;
	pcoord.x = (coord.x - (coord.x < 0)) / 2;
	pcoord.y = (coord.y - (coord.y < 0)) / 2;
	pcoord.z = (coord.z - (coord.z < 0)) / 2;

	return pcoord;
}

CellCoord child_coord(const CellCoord coord, uint32_t i)
{
	assert(i < 8);

	CellCoord ccoord;

	ccoord.lod = coord.lod - 1;
	ccoord.x = 2 * coord.x + ((i >> 0) & 1);
	ccoord.y = 2 * coord.y + ((i >> 1) & 1);
	ccoord.z = 2 * coord.z + ((i >> 2) & 1);

	return ccoord;
}

Mesh *MeshGrid::get_cell(CellCoord ccoord)
{
	uint32_t *p = cell_table.get(ccoord);
	return p ? &cells[*p] : NULL;
}

unsigned MeshGrid::get_children(CellCoord pcoord, Mesh *children[8])
{
	unsigned child_count = 0;

	for (int i = 0; i < 8; ++i) {
		CellCoord ccoord = child_coord(pcoord, i);
		uint32_t *p = cell_table.get(ccoord);
		if (p) {
			children[child_count] = &cells[*p];
			child_count++;
		}
	}

	return child_count;
}

MeshGrid::MeshGrid(Vec3 base, float step, uint32_t max_level, float err_tol)
    : base{base}, step{step}, levels{max_level + 1}, err_tol{err_tol},
      cell_offsets(levels), cell_counts(levels),
      cell_table(1 << (2 * levels + 4))
{
}

uint32_t MeshGrid::get_triangle_count(uint32_t level)
{
	if (level >= levels)
		return 0;

	uint32_t count = 0;

	for (uint32_t i = 0; i < cell_counts[level]; i++) {
		count += cells[i + cell_offsets[level]].index_count;
	}

	return (count / 3);
}

uint32_t MeshGrid::get_vertex_count(uint32_t level)
{
	if (level >= levels)
		return 0;

	uint32_t count = 0;

	for (uint32_t i = 0; i < cell_counts[level]; i++) {
		count += cells[i + cell_offsets[level]].vertex_count;
	}

	return (count);
}

static CellCoord block_base_coord(const CellCoord coord)
{
	CellCoord base;
	base.lod = coord.lod;
	base.x = ((coord.x + 1) & 0xFFFE) - 1;
	base.y = ((coord.y + 1) & 0xFFFE) - 1;
	base.z = ((coord.z + 1) & 0xFFFE) - 1;
	return base;
}

struct MeshGridBuilder {
	MeshGrid &mg;
	int num_threads;
	TArray<CellCoord> todo_blocks;
	int next_index = 0;
	pthread_mutex_t block_mutex;
	MeshGridBuilder(MeshGrid &mg, int num_thread);
	~MeshGridBuilder();
	void build_level(int level);
	void build_parent_cell(CellCoord pcoord);
	void build_block(CellCoord bcoord);
	bool build_next_block(void);
};

MeshGridBuilder::MeshGridBuilder(MeshGrid &mg, int num_threads)
    : mg{mg}, num_threads{num_threads}
{
	pthread_mutex_init(&block_mutex, NULL);
}

MeshGridBuilder::~MeshGridBuilder() { pthread_mutex_destroy(&block_mutex); }

void *build_blocks(void *args)
{
	struct MeshGridBuilder *builder = (struct MeshGridBuilder *)args;

	while (builder->build_next_block()) {
	}
	return NULL;
}

bool MeshGridBuilder::build_next_block()
{
	// MUTEX LOCK protection for next_index
	pthread_mutex_lock(&block_mutex);
	unsigned todo_index = next_index;
	if (todo_index >= todo_blocks.size) {
		pthread_mutex_unlock(&block_mutex);
		return false;
	}
	next_index += 1;
	pthread_mutex_unlock(&block_mutex);
	// MUTEX UNLOCK protection for next_index

	build_block(todo_blocks[todo_index]);

	return true;
}

void MeshGridBuilder::build_level(int level)
{

	/* Base level is built elsewhere */
	assert(level > 0);

	/* Pre allocate sufficently to avoid reallocation during of
	 * these during fct call, this allows to avoid some costly
	 * mutex lock.
	 */
	size_t max_cells = mg.cell_coords.size + mg.cell_counts[level - 1];
	mg.cell_coords.reserve(max_cells);
	mg.cells.reserve(max_cells);
	mg.cell_table.reserve(max_cells);

	/* Initiate, cell_counts[level] will be incremented during process */
	mg.cell_offsets[level] = mg.cells.size;
	mg.cell_counts[level] = 0;

	todo_blocks.reserve(mg.cell_counts[level - 1]);
	CellTable recorded_blocks(mg.cell_counts[level - 1]);

	/* Discover blocks and record them for later treatment */
	todo_blocks.clear();
	CellCoord *cell_coords = &mg.cell_coords[mg.cell_offsets[level - 1]];
	int num_cells = mg.cell_counts[level - 1];
	for (int i = 0; i < num_cells; ++i) {
		CellCoord ccoord = cell_coords[i];
		CellCoord pcoord = parent_coord(ccoord);
		CellCoord bcoord = block_base_coord(pcoord);
		if (!recorded_blocks.get(bcoord)) {
			recorded_blocks.set_at(bcoord, todo_blocks.size);
			todo_blocks.push_back(bcoord);
		}
	}

	/* Pre allocate an upper bound on the number of new indices and
	 * vertices in MBuf. This eases the parallel process of building
	 * blocks by avoiding a lot of otherwise necessary mutex locks.
	 */
	uint32_t alloc_idx = mg.get_triangle_count(level - 1) * 3;
	uint32_t alloc_vtx = mg.get_vertex_count(level - 1);
	mg.data.reserve_indices(mg.next_index_offset + alloc_idx);
	mg.data.reserve_vertices(mg.next_vertex_offset + alloc_vtx);

	/* Launch (and then join) threads to build separate blocks */
	pthread_t thread[num_threads];
	next_index = 0;
	for (int i = 0; i < num_threads; ++i) {
		pthread_create(&thread[i], NULL, build_blocks, (void *)this);
	}
	for (int i = 0; i < num_threads; ++i) {
		pthread_join(thread[i], NULL);
	}
}

void MeshGrid::build_from_mesh(const MBuf &src, const Mesh &mesh,
			       int num_threads)
{
	data.vtx_attr = src.vtx_attr | VtxAttr::MAP;

	init_from_mesh(src, mesh);
	printf("Number of cells at level 0 : %d\n", cell_counts[0]);

	MeshGridBuilder builder(*this, num_threads);

	for (uint32_t level = 1; level < levels; level++) {

		builder.build_level(level);
		printf("Number of cells at level %d : %d\n", level,
		       cell_counts[level]);
		printf("Number of triangles at level %d  : %d (ratio : %f)\n",
		       level, get_triangle_count(level),
		       (float)get_triangle_count(level) /
			   get_triangle_count(level - 1));
	}
}

void MeshGrid::init_from_mesh(const MBuf &src, const Mesh &mesh)
{
	cell_offsets[0] = 0;
	cell_counts[0] = 0;

	/* Loop over triangles and fill triangle to cell table. */

	size_t tri_count = mesh.index_count / 3;
	TArray<uint32_t> tri_idx_to_cell_idx(tri_count);

	CellCoord cell_coord{{0, 0, 0, 0}};
	float inv_step = 1.f / step;
	for (size_t tri_idx = 0; tri_idx < tri_count; ++tri_idx) {
		const uint32_t *indices = src.indices + mesh.index_offset;
		const Vec3 *positions = src.positions + mesh.vertex_offset;

		const Vec3 v1 = positions[indices[3 * tri_idx + 0]];
		const Vec3 v2 = positions[indices[3 * tri_idx + 1]];
		const Vec3 v3 = positions[indices[3 * tri_idx + 2]];

		const Vec3 bary = (v1 + v2 + v3) * (1.f / 3.f);

		point_to_cell_coord(cell_coord, bary, base, inv_step);

		uint32_t cell_idx;
		uint32_t *p = cell_table.get(cell_coord);
		if UNLIKELY (!p) {
			cell_idx = cell_table.size();
			cell_table.set_at(cell_coord, cell_idx);
			cells.push_back(Mesh{0, 0, 0, 0});
			cell_coords.push_back(cell_coord);
			cell_counts[0]++;
		} else {
			cell_idx = *p;
		}

		cells[cell_idx].index_count += 3;
		tri_idx_to_cell_idx[tri_idx] = cell_idx;
	}

	/* Compute index offsets and total index count */

	size_t total_index_count = 0;

	for (size_t cell_idx = 0; cell_idx < cells.size; ++cell_idx) {
		cells[cell_idx].index_offset = total_index_count;
		total_index_count += cells[cell_idx].index_count;
	}

	assert(total_index_count == mesh.index_count);

	/* Reserve indices */
	data.reserve_indices(total_index_count);

	/* Reset cells index_count to zero momentarily */
	for (size_t cell_idx = 0; cell_idx < cells.size; ++cell_idx) {
		cells[cell_idx].index_count = 0;
	}

	/* Reorder indices according to cells (and re-establish index counts) */
	for (size_t tri_idx = 0; tri_idx < tri_count; ++tri_idx) {
		Mesh &cell = cells[tri_idx_to_cell_idx[tri_idx]];

		size_t dst_offset = cell.index_offset + cell.index_count;
		size_t src_offset = mesh.index_offset + 3 * tri_idx;
		uint32_t *dst_idx = data.indices + dst_offset;
		uint32_t *src_idx = src.indices + src_offset;
		memcpy(dst_idx, src_idx, 3 * sizeof(uint32_t));

		cell.index_count += 3;
	}

	/* Clean-up */
	tri_idx_to_cell_idx.clear();

	/* Reserve vertices (total number is unknown at this point) */
	size_t vertex_count_guess = mesh.vertex_count + mesh.vertex_count / 4;
	data.reserve_vertices(vertex_count_guess);

	/**
	 * Compute max index_count in cells, this yields also an upper bound
	 * on max vertex_count per cell
	 */
	size_t max_index_count = 0;
	for (size_t cell_idx = 0; cell_idx < cells.size; ++cell_idx) {
		if (cells[cell_idx].index_count > max_index_count) {
			max_index_count = cells[cell_idx].index_count;
		}
	}

	HashTable<uint32_t, uint32_t> idx_remap(max_index_count);
	size_t total_vertex_count = 0;
	for (size_t cell_idx = 0; cell_idx < cells.size; ++cell_idx) {
		Mesh &cell = cells[cell_idx];
		cell.vertex_offset = total_vertex_count;

		data.reserve_vertices(total_vertex_count + cell.index_count);

		assert(cell.vertex_count == 0);

		uint32_t *cell_indices = data.indices + cell.index_offset;

		for (size_t i = 0; i < cell.index_count; ++i) {
			uint32_t old_idx = cell_indices[i];
			uint32_t new_idx;
			uint32_t *p = idx_remap.get(old_idx);
			if (!p) {
				new_idx = cell.vertex_count;
				idx_remap.set_at(old_idx, new_idx);
				size_t idx_s = mesh.vertex_offset + old_idx;
				size_t idx_d = total_vertex_count;
				copy_vertices(data, idx_d, src, idx_s, 1);
				cell.vertex_count++;
				total_vertex_count++;
			} else {
				new_idx = *p;
			}
			cell_indices[i] = new_idx;
		}
		/* Clear idx_remap for use with next cell */
		idx_remap.clear();
	}

	next_index_offset = total_index_count;
	next_vertex_offset = total_vertex_count;
}

void MeshGridBuilder::build_block(CellCoord bcoord)
{

	Mesh *children[8][8];

	/* For each parent cell in the block, the number of existing childs */
	uint32_t child_count[8];

	/* For each parent cell in the block, the number of indices and
	 * vertices prior to simplification */
	uint32_t idx_count[8];
	uint32_t vtx_count[8];

	/* For each parent cell in the block, its offsets in the index
	 * and vertex arrays. */
	uint32_t idx_offset[8];
	uint32_t vtx_offset[8];

	/* Some maximums for allocation bounds */
	uint32_t max_idx_count = 1;
	uint32_t max_vtx_count = 1;

	/* Count children for each parent cell in the block */
	for (uint32_t i = 0; i < 8; ++i) {
		CellCoord pcoord = bcoord;
		pcoord.x += (i >> 0) & 1;
		pcoord.y += (i >> 1) & 1;
		pcoord.z += (i >> 2) & 1;
		child_count[i] = mg.get_children(pcoord, children[i]);
	}

	/* Compute vtx and idx counts prior to simplification */
	uint32_t total_idx_count = 0;
	uint32_t total_vtx_count = 0;
	for (uint32_t i = 0; i < 8; i++) {
		idx_count[i] = 0;
		vtx_count[i] = 0;
		for (uint32_t j = 0; j < child_count[i]; ++j) {
			idx_count[i] += children[i][j]->index_count;
			vtx_count[i] += children[i][j]->vertex_count;
		}
		total_idx_count += idx_count[i];
		total_vtx_count += vtx_count[i];
		max_idx_count = MAX(max_idx_count, idx_count[i]);
		max_vtx_count = MAX(max_vtx_count, vtx_count[i]);
		idx_offset[i] =
		    (i == 0) ? 0 : idx_offset[i - 1] + idx_count[i - 1];
		vtx_offset[i] =
		    (i == 0) ? 0 : vtx_offset[i - 1] + vtx_count[i - 1];
	}

	/* Prepare tmp structures for simplification */

	MBuf blk_data;
	blk_data.vtx_attr = mg.data.vtx_attr;
	blk_data.reserve_indices(total_idx_count);
	blk_data.reserve_vertices(total_vtx_count + 1);

	VertexTable blk_table(total_vtx_count + 16, &blk_data,
			      blk_data.vtx_attr);

	TArray<uint32_t> blk_remap(total_vtx_count);
	for (uint32_t l = 0; l < total_vtx_count; ++l) {
		blk_remap[l] = ~0u;
	}

	Mesh blk_mesh{0, 0, 0, 0};

	/* Join meshes */

	uint32_t *remap = &blk_remap[0];
	for (uint32_t i = 0; i < 8; i++) {
		for (uint32_t j = 0; j < child_count[i]; ++j) {
			join_mesh_from_vertices(blk_mesh, blk_data,
						*children[i][j], mg.data,
						blk_table, remap);
			remap += children[i][j]->vertex_count;
		}
	}

	for (uint32_t k = 0; k < total_vtx_count; k++) {
		assert(blk_remap[k] < blk_mesh.vertex_count);
	}

	/* Simplify group */

	TArray<uint32_t> simp_remap(blk_mesh.vertex_count);
	TArray<uint32_t> trash(blk_mesh.index_count);
	float target_err = mg.err_tol * mg.step * (1 << bcoord.lod);
	meshopt_simplify_mod(
	    trash.data, simp_remap.data, blk_data.indices, blk_mesh.index_count,
	    (const float *)blk_data.positions, blk_mesh.vertex_count,
	    3 * sizeof(float), blk_mesh.index_count / 4, target_err, NULL);

	/* Update blk_mesh indices after simplification */

	for (uint32_t k = 0; k < blk_mesh.index_count; k++) {
		uint32_t *idx = blk_data.indices + blk_mesh.index_offset;
		assert(idx[k] < blk_mesh.vertex_count);
		idx[k] = simp_remap[idx[k]];
		assert(idx[k] < blk_mesh.vertex_count);
	}

	/* Compose simp_remap and blk_remap */

	for (uint32_t k = 0; k < total_vtx_count; k++) {
		blk_remap[k] = simp_remap[blk_remap[k]];
		assert(blk_remap[k] < blk_mesh.vertex_count);
	}

	/* Split each parent cell mesh from the block and copy it back
	 * to the mesh grid buffer data. */

	/* A second temp MBuf is allocated to spend less time inside mutex. */
	MBuf pdata;
	pdata.vtx_attr = mg.data.vtx_attr;
	pdata.reserve_indices(max_idx_count);
	pdata.reserve_vertices(max_vtx_count + 1);
	/* We recycle blk_table for pdata */
	blk_table.set_mesh_data(&pdata);

	TArray<uint32_t> split_remap(blk_mesh.vertex_count);

	for (uint32_t i = 0; i < 8; i++) {
		/* Go forward if this parent cell has no child */
		if (!child_count[i])
			continue;

		Mesh pmesh{0, 0, 0, 0};

		/* Clear blk_table and split remap */
		blk_table.clear();
		for (uint32_t l = 0; l < blk_mesh.vertex_count; ++l) {
			split_remap[l] = ~0u;
		}

		/* Perform the split : */

		/* 1) Copy relevant vertices and fill split_remap */
		pmesh.vertex_count = copy_unique_vertices(
		    pdata, 0, blk_data, &blk_remap[vtx_offset[i]], vtx_count[i],
		    blk_table, &split_remap[0]);

		/* If at max_level, write our own parent as ourself */
		if (bcoord.lod == (int16_t)(mg.levels - 1)) {
			for (uint32_t i = 0; i < pmesh.vertex_count; ++i) {
				pdata.remap[i] = i;
			}
		}

		/* 2) Copy relevant indices, skipping triangles that have
		 *    been collapsed by simplification */
		uint32_t p_idx_count = 0;
		uint32_t *tri_idx = blk_data.indices + idx_offset[i];
		for (uint32_t k = 0; k < idx_count[i]; k += 3) {
			uint32_t i0 = split_remap[tri_idx[k + 0]];
			uint32_t i1 = split_remap[tri_idx[k + 1]];
			uint32_t i2 = split_remap[tri_idx[k + 2]];
			assert(i0 != ~0u && i1 != ~0u && i2 != ~0u);
			if (i0 == i1 || i1 == i2 || i0 == i2)
				continue;
			pdata.indices[p_idx_count + 0] = i0;
			pdata.indices[p_idx_count + 1] = i1;
			pdata.indices[p_idx_count + 2] = i2;
			p_idx_count += 3;
		}
		pmesh.index_count = p_idx_count;

		/* 3) Record children cells parent map */
		uint32_t *src_idx = &blk_remap[vtx_offset[i]];
		for (uint32_t j = 0; j < child_count[i]; ++j) {
			const Mesh *cmesh = children[i][j];
			uint32_t *tgt_idx =
			    mg.data.remap + cmesh->vertex_offset;
			for (uint32_t k = 0; k < cmesh->vertex_count; k++) {
				tgt_idx[k] = split_remap[src_idx[k]];
				assert(tgt_idx[k] < pmesh.vertex_count);
			}
			src_idx += cmesh->vertex_count;
		}

		/* Write parent cell to mesh grid (protected by mutex) */
		CellCoord pcoord = bcoord;
		pcoord.x += (i >> 0) & 1;
		pcoord.y += (i >> 1) & 1;
		pcoord.z += (i >> 2) & 1;
		pthread_mutex_lock(&block_mutex);
		{
			pmesh.index_offset = mg.next_index_offset;
			pmesh.vertex_offset = mg.next_vertex_offset;
			uint32_t cell_idx = mg.cell_table.size();
			mg.cell_table.set_at(pcoord, cell_idx);
			mg.cells.push_back(pmesh);
			mg.cell_coords.push_back(pcoord);
			mg.cell_counts[pcoord.lod]++;
			mg.next_index_offset += pmesh.index_count;
			mg.next_vertex_offset += pmesh.vertex_count;
		}
		pthread_mutex_unlock(&block_mutex);

		copy_indices(mg.data, pmesh.index_offset, pdata, 0,
			     pmesh.index_count, 0);

		copy_vertices(mg.data, pmesh.vertex_offset, pdata, 0,
			      pmesh.vertex_count, 0);
	}
}

void MeshGridBuilder::build_parent_cell(CellCoord pcoord)
{
	Mesh *target_mesh = mg.get_cell(pcoord);
	assert(target_mesh);

	/* Get children */
	Mesh *children[8];
	unsigned child_count = mg.get_children(pcoord, children);
	assert(child_count);

	/* Compute vtx and idx counts prior to simplification */
	uint32_t idx_count = 0;
	uint32_t vtx_count = 0;
	for (uint32_t i = 0; i < child_count; ++i) {
		idx_count += children[i]->index_count;
		vtx_count += children[i]->vertex_count;
	}

	/* Prepare structures for simplification */

	MBuf tmp_data;
	tmp_data.vtx_attr = mg.data.vtx_attr;
	tmp_data.reserve_indices(idx_count);
	tmp_data.reserve_vertices(vtx_count + 1);

	VertexTable tmp_table(vtx_count + 16, &tmp_data, tmp_data.vtx_attr);

	TArray<uint32_t> remap(vtx_count);

	Mesh tmp_mesh{0, 0, 0, 0};

	/* Join meshes */
	uint32_t *tmp_remap = remap.data;
	for (uint32_t i = 0; i < child_count; ++i) {
		join_mesh_from_indices(tmp_mesh, tmp_data, *children[i],
				       mg.data, tmp_table, tmp_remap);
		tmp_remap += children[i]->vertex_count;
	}

	/* Simplify group */
	TArray<uint32_t> remap2(tmp_mesh.vertex_count);
	tmp_mesh.index_count = meshopt_simplify_mod(
	    tmp_data.indices, remap2.data, tmp_data.indices,
	    tmp_mesh.index_count, (const float *)tmp_data.positions,
	    tmp_mesh.vertex_count, 3 * sizeof(float), tmp_mesh.index_count / 4,
	    1, NULL);

	/* Update remap after simplification */
	for (uint32_t i = 0; i < vtx_count; i++) {
		remap[i] = remap2[remap[i]];
	}

	/* Write back group to mesh_grid */
	tmp_table.clear();
	tmp_table.set_mesh_data(&mg.data);
	target_mesh->index_offset = mg.next_index_offset;
	target_mesh->vertex_offset = mg.next_vertex_offset;
	mg.data.reserve_indices(mg.next_index_offset + idx_count);
	mg.data.reserve_vertices(mg.next_vertex_offset + vtx_count + 1);
	join_mesh_from_indices(*target_mesh, mg.data, tmp_mesh, tmp_data,
			       tmp_table, remap2.data);
	assert(idx_count >= target_mesh->index_count &&
	       vtx_count >= target_mesh->vertex_count);
	mg.next_index_offset += target_mesh->index_count;
	mg.next_vertex_offset += target_mesh->vertex_count;

	/* Update remap after write back */
	for (uint32_t i = 0; i < vtx_count; i++) {
		remap[i] = remap2[remap[i]];
	}

	/* Update child vertices parent map using remap */
	tmp_remap = remap.data;
	for (unsigned i = 0; i < child_count; i++) {
		uint32_t *p = mg.data.remap + children[i]->vertex_offset;
		for (uint32_t j = 0; j < children[i]->vertex_count; j++) {
			p[j] = tmp_remap[j];
		}
		tmp_remap += children[i]->vertex_count;
	}
}

void MeshGrid::select_cells_from_view_point(const Vec3 vp, float kappa,
					    const float *pvm,
					    TArray<uint32_t> &to_draw,
					    TArray<uint32_t> &parents)
{
	struct Candidate {
		uint32_t idx;
		uint32_t parent_idx;
		uint32_t check_visibility;
	};

	/* TODO avoid malloc in hot rendering loop */
	TArray<Candidate> to_visit;

	/* Load max level cell(s) */
	for (uint32_t i = 0; i < cell_counts[levels - 1]; i++) {
		uint32_t idx = cell_offsets[levels - 1] + i;
		Candidate candi = {idx, idx, 1};
		to_visit.push_back(candi);
	}

	size_t visited = 0;

	while (visited < to_visit.size) {
		Candidate candi = to_visit[visited++];
		CellCoord coord = cell_coords[candi.idx];

		/* Frustum cull */
		int visible = 1;
		if (pvm && candi.check_visibility) {
			Aabb bbox;
			bbox.min.x = base.x + coord.x * step * (1 << coord.lod);
			bbox.min.y = base.y + coord.y * step * (1 << coord.lod);
			bbox.min.z = base.z + coord.z * step * (1 << coord.lod);
			bbox.max.x = bbox.min.x + step * (1 << coord.lod);
			bbox.max.y = bbox.min.y + step * (1 << coord.lod);
			bbox.max.z = bbox.min.z + step * (1 << coord.lod);
			visible = is_visible(bbox, pvm);
			if (!visible)
				continue;
		}

		/* Distance based LOD */
		if (coord.lod == 0 || cell_view_ratio(vp, coord) > kappa) {
			to_draw.push_back(candi.idx);
			parents.push_back(candi.parent_idx);
		} else {
			uint32_t check_vis = visible != 2;
			for (int i = 0; i < 8; ++i) {
				CellCoord ccoord = child_coord(coord, i);
				uint32_t *p = cell_table.get(ccoord);
				if (p) {
					to_visit.push_back(
					    {*p, candi.idx, check_vis});
				}
			}
		}
	}
}

float MeshGrid::cell_view_ratio(Vec3 vp, CellCoord coord)
{
	Vec3 diff = (vp - base) * (1.f / (step * (1 << coord.lod)));

	float rat = 0;

	rat = coord.x - diff.x > rat ? coord.x - diff.x : rat;
	rat = diff.x - (coord.x + 1) > rat ? diff.x - (coord.x + 1) : rat;
	rat = coord.y - diff.y > rat ? coord.y - diff.y : rat;
	rat = diff.y - (coord.y + 1) > rat ? diff.y - (coord.y + 1) : rat;
	rat = coord.z - diff.z > rat ? coord.z - diff.z : rat;
	rat = diff.z - (coord.z + 1) > rat ? diff.z - (coord.z + 1) : rat;

	return rat;
}
