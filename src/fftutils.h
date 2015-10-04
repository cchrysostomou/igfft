#ifndef FFT_H
#define FFT_H

#include<cmath>
#include<complex>
#include<iostream>
#include<vector>

#define Vec(a, b) std::vector<__typeof(*(a))> ((a), (a)+(b))

// allow easy change to float or long double
#define USE_DOUBLE

#ifdef USE_FLOAT
typedef float complex complex_t;
typedef float real_t;
#define cexp cexpf
#define exp expf
#endif


#ifdef USE_DOUBLE
typedef complex<double> complex_t;
typedef double real_t;
#endif



#endif