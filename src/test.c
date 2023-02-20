#include <stdio.h>
#include <stdlib.h>

#include "arc_raster.h"

int main(int argc, char **argv)
{
	struct ArcRaster raster;
	raster.data = (float *)malloc(1000 * 1000 * sizeof(float));

	if (read_ArcRaster(argv[1], &raster)) return -1;
	
	double x = atof(argv[2]);
	double y = atof(argv[3]);

	printf("Altitude : %f\n", ArcRaster_elev(x, y, &raster));

	return 0;
}
