#ifndef TXT_IO_FAST_H
#define TXT_IO_FAST_H

#ifdef __cplusplus
extern "C" {
#endif

unsigned get_next_unsigned(const char **str);

int get_next_int(const char **str);

double get_next_real(const char **str);

void get_next_string(const char **str, int n, char *out);

void jump_next(const char **str, int n);

#ifdef __cplusplus
}
#endif

#endif
