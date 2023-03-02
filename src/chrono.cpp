#include <sys/time.h>
#include <stdio.h>

#include <chrono.h>

struct timeval tv0, tv1;

void timer_start()
{
	gettimeofday(&tv0, NULL);
}

unsigned int timer_stop(const char *str)
{
	gettimeofday(&tv1, NULL);
	unsigned int mus = 1000000 * (tv1.tv_sec - tv0.tv_sec);
	mus += (tv1.tv_usec - tv0.tv_usec);
	if (str[0])
	{
		printf("Timer %s: ", str);
		if  (mus >= 1000000) 
		{
			printf("%.3f s\n", (float)mus / 1000000);
		}
		else 
		{
			printf("%.3f ms\n", (float)mus / 1000);
		}
	}
	return (mus);
}

