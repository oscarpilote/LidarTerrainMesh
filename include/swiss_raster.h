#pragma once

int read_swiss_geotiff(float *dest, int *pix_x, int *pix_y, const char *fname,
		       bool size_only = false);

