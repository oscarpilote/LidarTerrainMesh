struct ArcRaster {
	unsigned ncols;
	unsigned nrows;
	double xllcenter;
	double yllcenter;
	double cellsize;
	double NODATA_value;
	float *data;
};

int read_ArcRaster(const char *filename, struct ArcRaster* r);

/* Bilinear interpolation + closest match for out of bounds 
 * WARN: Do not take care of NODATA_value.
 */
float ArcRaster_elev(double x, double y, const struct ArcRaster *r); 

