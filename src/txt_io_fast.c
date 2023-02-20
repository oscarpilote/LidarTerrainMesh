#include "txt_io_fast.h"

static inline const char* skip_separators(const char *str)
{
	while(*str <= ' ') {str++;}
	return str;
}

unsigned get_next_unsigned(const char **str) {
	unsigned ret=0;
	const char *s = *str;
	s = skip_separators(s);
	while (*s >= '0') {
		ret*=10;
		ret+=(*s-'0');
		++s;
	}
	s = skip_separators(s);
	*str = s;
	return ret;
}

int get_next_int(const char **str) {
	int sign = 1;
	unsigned ret=0;
	const char *s = *str;
	s = skip_separators(s);
	if (*s == '-')
	{
		sign = -1;
		++s;
	}
	while (*s >= '0')
	{
		ret*=10;
		ret+=(*s-'0');
		++s;
	}
	s = skip_separators(s);
	*str = s;
	return sign * ret;
}

double get_next_real(const char **str) 
{
	int sign = 1;
	double ret = 0;
	const char *s = *str;
	s = skip_separators(s);
	if (*s == '-') {
		sign = -1;
		++s;
	}
	while (*s >= '0') {
		ret *= 10;
		ret += (*s - '0');
		++s;
	}
	if (*s == '.')
	{
		++s;
		double deci = 1;
		while (*s >= '0') 
		{
			ret *= 10;
			ret += (*s - '0');
			deci *= 10;
			++s;
		}
		ret /= deci;
	}
	s = skip_separators(s);
	*str = s;
	return sign * ret;
}

void get_next_string(const char **str, int n, char *out)
{
	const char *s = *str;
	int len = n;
	s = skip_separators(s);
	while (*s > ' ' && len--) {
		*out++ = *s++;
	}
	*out = '\0';
	s = skip_separators(s);
	*str = s;
}

void jump_next(const char **str, int n)
{
	const char *s = *str;
	while (n--) {
		while (*s <= ' ') {++s;}
		while (*s > ' ') {++s;}
	}
	while (*s <= ' ') {++s;}
	*str = s;
}
