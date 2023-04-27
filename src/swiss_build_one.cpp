#include <stdlib.h>
#include <string.h>

#include "multigrid.h"
#include "swiss_mgrid.h"

int main(int argc, char **argv)
{
	char dummy[5] = {0};
	strncpy(dummy, argv[1], 4);
	int nx = atoi(dummy);
	strncpy(dummy, argv[1] + 5, 4);
	int ny = atoi(dummy);
	Pyramid p = {{nx, ny, 0}, 2};
	init_mgrid(argv[1], argv[2], p);
	return (0);
}
