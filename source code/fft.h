#ifndef FFT_H
#define FFT_H

#include <vector>
#include <cmath>
#include <cstdio>
//#include <complex.h>
#include "fftw3.h"


#define Vec(a, b) std::vector<__typeof(*(a))> ((a), (a)+(b))

// allow easy change to float or long double
//#define USE_FLOAT
//#define USE_DOUBLE
#define USE_FFTWCOMPLEX

#ifdef USE_FLOAT
typedef float complex_t;
typedef float real_t;
#define cexp cexpf
#define exp expf
#endif

#ifdef USE_DOUBLE
typedef double _Complex complex_t;
typedef double real_t;
#endif

#ifdef USE_FFTWCOMPLEX
typedef fftw_complex complex_t;
typedef double real_t;
#endif

//#define DEBUG

#endif
