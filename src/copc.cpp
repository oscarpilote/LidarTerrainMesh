#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "lazperf/readers.hpp"

#include "aabb.h"
#include "array.h"
#include "hash.h"
#include "hash_table.h"
#include "math_utils.h"
//#include "las_read.h"

#include "copc.h"

#define LAS_HEADER_SIZE 375
#define VLR_HEADER_SIZE 54

struct CopcInfo {
	double center_x;
	double center_y;
	double center_z;
	double halfsize;
	double spacing;
	uint64_t root_hier_offset;
	uint64_t root_hier_size;
	double gpstime_minimum;
	double gpstime_maximum;
	uint64_t reserved[11];
};

struct CellInfo {
	uint64_t offset;
	uint32_t byte_size;
	int32_t point_count;
};

union CellKey {
	struct {
		int32_t level;
		int32_t x;
		int32_t y;
		int32_t z;
	};
	struct {
		uint64_t low;
		uint64_t high;
	};
};

static CellKey child(CellKey &key, int i)
{
	CellKey child;
	child.level = key.level + 1;
	child.x = 2 * key.x + ((i & 1) >> 0);
	child.y = 2 * key.y + ((i & 2) >> 1);
	child.z = 2 * key.z + ((i & 4) >> 2);
	return child;
}

struct CellKeyHasher {
	static constexpr CellKey empty_key = {{-1, 0, 0, 0}};
	size_t hash(CellKey key) const;
	bool is_empty(CellKey key) const;
	bool is_equal(CellKey key1, CellKey key2) const;
};

inline size_t CellKeyHasher::hash(CellKey key) const
{
	size_t hash = murmur2_64(0, key.low);
	return murmur2_64(hash, key.high);
}

inline bool CellKeyHasher::is_empty(CellKey key) const
{
	return key.level == -1;
}

inline bool CellKeyHasher::is_equal(CellKey key1, CellKey key2) const
{
	return key1.low == key2.low && key1.high == key2.high;
}

struct VLRHeader {
	uint16_t reserved;
	char user_id[16];
	uint16_t record_id;
	uint16_t record_length_after_header;
	char description[32];
};

struct CopcReader {
	int point_type;
	int point_size;
	TVec3<double> cube_base;
	double cube_size;
	double spacing;
	FILE *f;
	HashTable<CellKey, CellInfo, CellKeyHasher> cells;
	TArray<CellInfo *> bbox_cells;
	TArray<char> scratch;
	uint32_t max_scratch_need;
	/* Methods */
	int open(const char *filename);
	void read_cells_info(size_t offset, uint32_t size);
	uint32_t set_target_bbox(const TAabb<double> &box);
	void read_cell(CellInfo *info, char *out);

	uint32_t bound_overlapping_count(const TAabb<double> &box,
					 CellKey key = {{0, 0, 0, 0}});
	uint32_t read_overlapping_cells(const TAabb<double> &box, char *out,
					CellKey key = {{0, 0, 0, 0}});
};

CopcReader *copc_init(const char *filename)
{
	CopcReader *copc = new CopcReader;
	if (copc->open(filename)) {
		delete copc;
		return nullptr;
	}
	return copc;
}

void copc_fini(CopcReader *copc) { delete copc; }

uint32_t copc_set_target_bbox(CopcReader *copc, const TAabb<double> &box)
{
	return copc->set_target_bbox(box);
}

int copc_cell_point_count(CopcReader *copc, uint32_t cell_idx)
{
	if (cell_idx < copc->bbox_cells.size)
		return copc->bbox_cells[cell_idx]->point_count;
	return 0;
}

int copc_read_cell(CopcReader *copc, uint32_t cell_idx, char *dst)
{
	if (cell_idx >= copc->bbox_cells.size)
		return -1;
	assert(dst);
	copc->read_cell(copc->bbox_cells[cell_idx], dst);
	return 0;
}

uint32_t copc_bound_inside(CopcReader *copc, const TAabb<double> &box)
{
	return copc->bound_overlapping_count(box);
}

uint32_t copc_load_inside(CopcReader *copc, const TAabb<double> &box, char *buf)
{
	return copc->read_overlapping_cells(box, buf);
}

int CopcReader::open(const char *filename)
{
	f = fopen(filename, "rb");
	if (!f)
		return -1;

	char buf[589];
	if (fread(buf, 549, 1, f) != 1) {
		fclose(f);
		f = nullptr;
		return -1;
	}

	if (strncmp(buf, "LASF", 4) != 0) {
		fclose(f);
		f = nullptr;
		printf("Not a LASF file\n");
		return -1;
	}

	struct VLRHeader *h = (struct VLRHeader *)(buf + LAS_HEADER_SIZE);
	if (strncmp(h->user_id, "copc", 4) != 0 || h->record_id != 1) {
		printf("Wrong COPC header\n");
		fclose(f);
		f = nullptr;
		return -1;
	}

	point_type = *(unsigned char *)(buf + 104) & 0x3F;
	point_size = *(unsigned char *)(buf + 105);

	struct CopcInfo *i =
	    (struct CopcInfo *)(buf + LAS_HEADER_SIZE + VLR_HEADER_SIZE);
	cube_base.x = i->center_x - i->halfsize;
	cube_base.y = i->center_y - i->halfsize;
	cube_base.z = i->center_z - i->halfsize;
	cube_size = 2 * i->halfsize;
	spacing = i->spacing;

	uint64_t root_offset = i->root_hier_offset;
	uint64_t root_size = (uint32_t)i->root_hier_size;
	max_scratch_need = 0;
	read_cells_info(root_offset, root_size);
	return 0;
}

void CopcReader::read_cells_info(size_t offset, uint32_t size)
{
	assert(f);
	assert(size % 32 == 0);
	int entries = size / 32;
	CellKey key;
	CellInfo info;
	for (int i = 0; i < entries; ++i) {
		fseek(f, offset, SEEK_SET);
		fread(&key, 16, 1, f);
		fread(&info, 16, 1, f);
		if (info.point_count < 0) {
			read_cells_info(info.offset, info.byte_size);
		} else {
			assert(!cells.get(key));
			cells.set_at(key, info);
			max_scratch_need =
			    MAX(max_scratch_need, info.byte_size);
		}
		offset += 32;
	}
}

uint32_t CopcReader::bound_overlapping_count(const TAabb<double> &box,
					     CellKey key)
{
	CellInfo *info = cells.get(key);
	if (!info)
		return 0;
	assert(info->point_count >= 0);

	TAabb<double> cbox;
	double scale = cube_size / (1 << key.level);
	cbox.min.x = cube_base.x + key.x * scale;
	cbox.min.y = cube_base.y + key.y * scale;
	cbox.min.z = cube_base.z + key.z * scale;
	cbox.max.x = cbox.min.x + scale;
	cbox.max.y = cbox.min.y + scale;
	cbox.max.z = cbox.min.z + scale;
	cbox &= box;
	if (cbox.is_empty()) {
		return 0;
	}

	uint32_t bound = info->point_count;
	for (int i = 0; i < 8; i++) {
		CellKey child_key = child(key, i);
		bound += bound_overlapping_count(box, child_key);
	}
	return bound;
}

void CopcReader::read_cell(CellInfo *info, char *out)
{
	scratch.reserve(max_scratch_need);
	fseek(f, info->offset, SEEK_SET);
	char *src = &scratch[0];
	fread(src, info->byte_size, 1, f);
	lazperf::reader::chunk_decompressor dec(point_type, 0, src);
	for (int i = 0; i < info->point_count; ++i) {
		dec.decompress(out);
		out += point_size;
	}
}

uint32_t CopcReader::set_target_bbox(const TAabb<double> &box)
{
	bbox_cells.clear();
	TArray<CellKey> to_visit;
	to_visit.push_back({{0, 0, 0, 0}});
	uint32_t visited = 0;
	while (visited < to_visit.size) {
		CellKey key = to_visit[visited];
		visited++;

		CellInfo *info = cells.get(key);
		/* If cell not present continue */
		if (!info)
			continue;

		/* If cell bbox does not intersect target bbxo continue */
		TAabb<double> cbox;
		double scale = cube_size / (1 << key.level);
		cbox.min.x = cube_base.x + key.x * scale;
		cbox.min.y = cube_base.y + key.y * scale;
		cbox.min.z = cube_base.z + key.z * scale;
		cbox.max.x = cbox.min.x + scale;
		cbox.max.y = cbox.min.y + scale;
		cbox.max.z = cbox.min.z + scale;
		// printf("%lf %lf %lf %lf %lf %lf\n", box.min.x, box.max.x,
		//        box.min.y, box.max.y, box.min.z, box.max.z);
		// printf("%lf %lf %lf %lf %lf %lf\n", cbox.min.x, cbox.max.x,
		//        cbox.min.y, cbox.max.y, cbox.min.z, cbox.max.z);
		cbox &= box;
		// printf("%lf %lf %lf %lf %lf %lf\n", cbox.min.x, cbox.max.x,
		//        cbox.min.y, cbox.max.y, cbox.min.z, cbox.max.z);
		if (cbox.is_empty())
			continue;

		/* Push into bbox_cells and push children for visiting */
		bbox_cells.push_back(info);
		for (int i = 0; i < 8; ++i) {
			CellKey child_key = child(key, i);
			to_visit.push_back(child_key);
		}
	}
	return bbox_cells.size;
}

uint32_t CopcReader::read_overlapping_cells(const TAabb<double> &box, char *out,
					    CellKey key)
{
	CellInfo *info = cells.get(key);
	if (!info)
		return 0;

	assert(info->point_count >= 0);

	TAabb<double> cbox;
	double scale = cube_size / (1 << key.level);
	cbox.min.x = cube_base.x + key.x * scale;
	cbox.min.y = cube_base.y + key.y * scale;
	cbox.min.z = cube_base.z + key.z * scale;
	cbox.max.x = cbox.min.x + scale;
	cbox.max.y = cbox.min.y + scale;
	cbox.max.z = cbox.min.z + scale;
	cbox &= box;
	if (cbox.is_empty()) {
		return 0;
	}

	uint32_t count = info->point_count;
	if (count > 0) {
		read_cell(info, out);
	}
	for (int i = 0; i < 8; i++) {
		CellKey child_key = child(key, i);
		count += read_overlapping_cells(box, out + count * point_size,
						child_key);
	}
	return count;
}

/*int main(int argc, char **argv)
{
	LasFileInfo info;
	las_read_info(argv[1], info);
	las_print_info(info);
	struct CopcReader copc;
	if (copc.open(argv[1]))
		return -1;
	TAabb<double> box;
	box.min = copc.cube_base;
	box.min.x += copc.cube_size * atof(argv[2]);
	box.max = copc.cube_base + copc.cube_size * TVec3<double>{1., 1., 1.};
	printf("%f %f %f %f\n", box.min.x, box.max.x, box.min.y, box.max.y);
	int bound = copc.bound_overlapping_count(box);
	char *out = (char *)malloc(bound * copc.point_size);
	copc.read_overlapping_cells(box, out);
	printf("Number of points read : %.2fM\n", (float)bound * 1e-6);
	return 0;
}*/
