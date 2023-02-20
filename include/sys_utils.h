#pragma once

/* C++ does not allow implicitit cast of void*   */
/* This will work with gcc or clang              */
#if defined(__cplusplus)
	#define VOIDSTARCAST(x) (typeof (x))
#else
	#define VOIDSTARCAST(x) 
#endif

/* Branch prediction */
#if defined(__GNUC__) || defined(__clang__)
	#define LIKELY(x) (__builtin_expect((x), 1))
	#define UNLIKELY(x) (__builtin_expect((x), 1))
	#define ASSUME(cond) do { if (!(cond)) __builtin_unreachable(); } while (0)
	#define ASSUME_ALIGNED(var, size) do { var = VOIDSTARCAST(var) __builtin_assume_aligned(var, size); } while (0)
#else
	#define LIKELY(x) (x)
	#define UNLIKELY(x) (x)
	#define ASSUME(cond)
	#define ASSUME_ALIGNED(var, n)
#endif

#define MALLOC_NUM(ptr, num) do {\
	ptr = VOIDSTARCAST(ptr)malloc((num) * sizeof (*ptr));\
	if (UNLIKELY(num != 0 && ptr == NULL)) abort();\
	} while(0)

#define MALLOC_SIZ(ptr, size) do {\
	ptr = VOIDSTARCAST(ptr)malloc((size));\
	if (UNLIKELY(size != 0 && ptr == NULL)) abort();\
	} while(0)

#define CALLOC_NUM(ptr, num) do {\
	ptr = VOIDSTARCAST(ptr)calloc((num), sizeof (*ptr));\
	if (UNLIKELY(num != 0 && ptr == NULL)) abort();\
	} while(0)

#define CALLOC_SIZ(ptr, size) do {\
	ptr = VOIDSTARCAST(ptr)calloc((size) / sizeof (*ptr), sizeof (*ptr));\
	if (UNLIKELY(size != 0 && ptr == NULL)) abort();\
	} while(0)

#define REALLOC_NUM(ptr, num) do {\
	ptr = VOIDSTARCAST(ptr)realloc(ptr, (num) * sizeof (*ptr));\
	if (UNLIKELY(num != 0 && ptr == NULL)) abort();\
	} while(0)

#define REALLOC_SIZ(ptr, size) do {\
	ptr = VOIDSTARCAST(ptr)realloc(ptr, (size));\
	if (UNLIKELY(size != 0 && ptr == NULL)) abort();\
	} while(0)

#define MEMFREE(ptr) do {free(ptr); ptr = NULL;} while (0)

// Remove old 16bit Windows.h macros near and far
#ifdef far
#undef far
#endif
#ifdef near
#undef near
#endif
