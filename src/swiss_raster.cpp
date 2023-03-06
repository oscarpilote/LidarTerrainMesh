#include "gdal/gdal_priv.h"
#include "stdio.h"

#include "swiss_raster.h"

int read_swiss_geotiff(float *dest, int *pix_x, int *pix_y, const char *fname,
		       bool size_only)
{
	{
		FILE *f;
		if ((f = fopen(fname, "rb")) == NULL)
			return -1;
		fclose(f);
	}

	GDALDatasetH hDataset;
	GDALAllRegister();
	hDataset = GDALOpen(fname, GA_ReadOnly);
	if (!hDataset)
		return (-1);
	GDALRasterBandH hBand;
	hBand = GDALGetRasterBand(hDataset, 1);
	int nx = GDALGetRasterBandXSize(hBand);
	int ny = GDALGetRasterBandYSize(hBand);
	if (!size_only) {
		for (int i = 0; i < ny; ++i) {
			float *target = dest + i * nx;
			int dummy =
			    GDALRasterIO(hBand, GF_Read, 0, i, nx, 1, target,
					 nx, 1, GDT_Float32, 0, 0);
			(void)dummy;
		}
	}

	*pix_x = nx;
	*pix_y = ny;

	GDALClose(hDataset);

	return (nx * ny != 0 ? 0 : -1);
}

