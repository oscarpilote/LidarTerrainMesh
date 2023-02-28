#pragma once

#define MAX_POWER 20

static const double POWER_10_POS[MAX_POWER] = {
    1.0e0,  1.0e1,  1.0e2,  1.0e3,  1.0e4,  1.0e5,  1.0e6,
    1.0e7,  1.0e8,  1.0e9,  1.0e10, 1.0e11, 1.0e12, 1.0e13,
    1.0e14, 1.0e15, 1.0e16, 1.0e17, 1.0e18, 1.0e19,
};

static const double POWER_10_NEG[MAX_POWER] = {
    1.0e0,   1.0e-1,  1.0e-2,  1.0e-3,	1.0e-4,	 1.0e-5,  1.0e-6,
    1.0e-7,  1.0e-8,  1.0e-9,  1.0e-10, 1.0e-11, 1.0e-12, 1.0e-13,
    1.0e-14, 1.0e-15, 1.0e-16, 1.0e-17, 1.0e-18, 1.0e-19,
};

static inline int is_whitespace(char c)
{
	return (c == ' ' || c == '\t' || c == '\r');
}

static inline int is_newline(char c) { return (c == '\n'); }

static inline int is_digit(char c) { return (c >= '0' && c <= '9'); }

static inline int is_exponent(char c) { return (c == 'e' || c == 'E'); }

static inline const char *skip_whitespace(const char *ptr)
{
	while (is_whitespace(*ptr))
		ptr++;

	return ptr;
}

static inline const char *skip_line(const char *ptr)
{
	while (!is_newline(*ptr++))
		;

	return ptr;
}

static const char *parse_int(const char *ptr, int *val)
{
	int sign;
	int num;

	ptr = skip_whitespace(ptr);

	if (*ptr == '-') {
		sign = -1;
		ptr++;
	} else {
		sign = +1;
	}

	num = 0;
	while (is_digit(*ptr))
		num = 10 * num + (*ptr++ - '0');

	*val = sign * num;

	return ptr;
}

static const char *parse_float(const char *ptr, float *val)
{
	double sign;
	double num;
	double fra;
	double div;
	int eval;
	const double *powers;

	ptr = skip_whitespace(ptr);

	switch (*ptr) {
	case '+':
		sign = 1.0;
		ptr++;
		break;

	case '-':
		sign = -1.0;
		ptr++;
		break;

	default:
		sign = 1.0;
		break;
	}

	num = 0.0;
	while (is_digit(*ptr))
		num = 10.0 * num + (double)(*ptr++ - '0');

	if (*ptr == '.')
		ptr++;

	fra = 0.0;
	div = 1.0;

	while (is_digit(*ptr)) {
		fra = 10.0 * fra + (double)(*ptr++ - '0');
		div *= 10.0;
	}

	num += fra / div;

	if (is_exponent(*ptr)) {
		ptr++;

		switch (*ptr) {
		case '+':
			powers = POWER_10_POS;
			ptr++;
			break;

		case '-':
			powers = POWER_10_NEG;
			ptr++;
			break;

		default:
			powers = POWER_10_POS;
			break;
		}

		eval = 0;
		while (is_digit(*ptr))
			eval = 10 * eval + (*ptr++ - '0');

		num *= (eval >= MAX_POWER) ? 0.0 : powers[eval];
	}

	*val = (float)(sign * num);

	return ptr;
}
