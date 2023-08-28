#pragma once

#include <sys/time.h>

struct Timer {
	bool launched = false;
	timeval tv0;
	timeval tv1;
	void start();
	unsigned int stop(const char *str = "");
};
