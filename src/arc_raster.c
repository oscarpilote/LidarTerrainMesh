#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include <math.h>

#include "txt_io_fast.h"
#include "arc_raster.h"

#define BUFFER_SIZE 65536

struct ArcRaster_reader {
	struct ArcRaster *raster;
	bool header_read;
	size_t num_data_read;
};

static void parse_buffer(const char *buffer, const char *end,
		struct ArcRaster_reader *reader)
{
	struct ArcRaster* r = reader->raster;

	const char *token = buffer;
	char tmp_str[256];

	if (!reader->header_read)
	{
		get_next_string(&token, 255, tmp_str);
		assert(strncmp(tmp_str, "ncols", 5) == 0);  
		r->ncols = get_next_unsigned(&token);
		get_next_string(&token, 255, tmp_str);
		assert(strncmp(tmp_str, "nrows", 5) == 0);  
		r->nrows = get_next_unsigned(&token);
		get_next_string(&token, 255, tmp_str);
		bool center = true;
		if (strncmp(tmp_str, "xllcenter", 9) != 0)
		{
			center = false;
			assert(strncmp(tmp_str, "xllcorner", 9) == 0);
		}
		r->xllcenter = get_next_real(&token);
		get_next_string(&token, 255, tmp_str);
		r->yllcenter = get_next_real(&token);
		get_next_string(&token, 255, tmp_str);
		r->cellsize = get_next_real(&token);
		if (!center)
		{
			r->xllcenter += r->cellsize / 2.0;
			r->yllcenter += r->cellsize / 2.0;
		}
		get_next_string(&token, 255, tmp_str);
		r->NODATA_value = get_next_real(&token);

		reader->header_read = true;
	}

	size_t num = reader->num_data_read; 
	while (token < end)
	{
		r->data[num++] = get_next_real(&token);
	}
	reader->num_data_read = num;
}

int read_ArcRaster(const char *filename, struct ArcRaster* raster)
{

	/* Open file */
	FILE *file = fopen(filename, "rb");
	if (!file)
		return -1;

	 /* Create buffer for reading file */
	char *buffer = (char*)malloc(2 * BUFFER_SIZE * sizeof(char));
	if (!buffer)
		return -1;

	struct ArcRaster_reader reader = {raster, false, 0};

	char *start = buffer;
	for (;;)
	{
		/* Read another buffer's worth from file */
		size_t read = fread(start, 1, BUFFER_SIZE, file);
		if (read == 0 && start == buffer)
			break;

		/* Ensure buffer ends in a newline */
		if (read < BUFFER_SIZE)
		{
			if (read == 0 || start[read - 1] != '\n')
				start[read++] = '\n';
		}
		char *end = start + read;


		/* Find last whitespace */
		char *last = end;
		while (last > buffer)
		{
			last--;
			if (*last == ' ' || *last == '\t' || *last == '\n')
				break;
		}


		/* Check there actually is a whitespace */
		if (*last != ' ' && *last != '\t' && *last != '\n')
			break;
		last++;


		/* Parse buffer */
		parse_buffer(buffer, last, &reader);


		/* Copy overflow for next buffer */
		size_t bytes = (size_t)(end - last);
		memmove(buffer, last, bytes);
		start = buffer + bytes;
	}


    /* Clean up */
    free(buffer);
    fclose(file);

    return 0;
}

float ArcRaster_elev(double x, double y, const struct ArcRaster* r)
{
	double cx = (x - r->xllcenter) / r->cellsize;
	int cl, cr;
	cl = floor(cx);
	cr = cl + 1;
	double alpha = cx - cl;
	/* Out of bounds check */
	if (cl < 0) 
	{
		cl = cr = 0;
	}
	else if (cr >= (int)r->ncols)
	{
		cl = cr = r->ncols - 1;
	}

	double ry = (y - r->yllcenter) / r->cellsize;
	int rb, rt;
	rb = floor(ry);
	rt = rb + 1;
	double beta = ry - rb;
	/* Out of bounds check */
	if (rb < 0) 
	{
		rb = rt = 0;
	}
	else if (rt >= (int)r->nrows)
	{
		rb = rt = r->nrows - 1;
	}
	/* first row is top */
	rb = r->nrows - 1 - rb;
	rt = r->nrows - 1 - rt;

	const float *arr = r->data;
	float res = (1 - alpha) * (1 - beta) * arr[rb * r->ncols + cl]
		  + (1 - alpha) * beta * arr[rt * r->ncols + cl]
		  + alpha * (1 - beta) * arr[rb * r->ncols + cr]
		  + alpha * beta * arr[rt * r->ncols + cr];

	return res;
}

