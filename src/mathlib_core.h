#ifndef ML_CORE_H
#define ML_CORE_H

#ifndef INTEGRAL_ITERATIONS
  #define INTEGRAL_ITERATIONS 1000000
#endif // INTEGRAL_ITERATIONS

#ifndef SIN_ITERATIONS
  #define SIN_ITERATIONS 4
#endif // SIN_ITERATIONS

#ifndef LN_ITERATIONS
  #define LN_ITERATIONS 1000
#endif // LN_ITERATIONS

#ifndef EXP_ITERATIONS
  #define EXP_ITERATIONS 50
#endif // EXP_ITERATIONS

#ifndef ARCSIN_ITERATIONS
  #define ARCSIN_ITERATIONS 33
#endif // ARCSIN_ITERATIONS

#ifndef ARCTAN_ITERATIONS
  #define ARCTAN_ITERATIONS 500
#endif // ARCTAN_ITERATIONS

// Included ONLY for malloc(), realloc() and free().
#include <stdlib.h>

// Taken from math.h -> Used to prevent any overhead from library such as trig functions.
/////////////////////////////////////////////////////////////////////////
#ifndef _HUGE_ENUF
    #define _HUGE_ENUF  1e+300  // _HUGE_ENUF*_HUGE_ENUF must overflow
#endif

#define INFINITY   ((float)(_HUGE_ENUF * _HUGE_ENUF))
#define HUGE_VAL   ((double)INFINITY)
#define HUGE_VALF  ((float)INFINITY)
#define HUGE_VALL  ((long double)INFINITY)
#ifndef _UCRT_NEGATIVE_NAN
// This operation creates a negative NAN adding a - to make it positive
#define NAN        (-(float)(INFINITY * 0.0F))
#else
// Keep this for backwards compatibility
#define NAN        ((float)(INFINITY * 0.0F))
#endif
//////////////////////////////////////////////////////////////////////////

#define MLDEF static inline
#define SMALL 0.000001f
#define ISINF(x) (x == INFINITY || x == -INFINITY)? 1 : 0
#define ISNAN(x) (x != x)? 1 : 0
#define CALL(f, x) (*f)(x)

typedef long double real;
typedef long int integer;
typedef unsigned long int natural;

typedef real (*function)(real n);

// Used for natural exponents. (positive integers or 0 ONLY.)
// For other real exponents, use pow_exp.
MLDEF real power(real x, natural n) {
  real result = 1;
  for(natural i = 0; i < n; ++i) {
    result *= x;
  }

  return result;
}

const real pi = 22.0 / 7.0;
const real e = 2.7182818284590452354;

// Uses recursion to calculate the factorial of a natural number.
MLDEF natural fact(natural n) {
  if(n == 0) return 1;

  return n * fact(n-1);
}

// Uses discrete_product to calculate the factorial of a natural number.
MLDEF natural fact_norecurse(natural n) {
  natural result = 1;

  for(int x = 1; x < n+1; ++x) {
    result *= x;
  }

  return result;
}

// Function that simply returns the input.
MLDEF real constant(real x) {return x;}

MLDEF real ml_abs(real x) {return (x >= 0)? x : (x*-1);}

MLDEF real ml_floor(real x) { return (real)((natural)x); }
MLDEF real ml_ceil(real x) { 
  real fractionalPart = x - ml_floor(x);

  return (fractionalPart == 0)? x : (ml_floor(x)+1);
}

// Natural logarithm -- logarithm with base e
MLDEF real ml_ln(real x) {
  real result = 0;

  for(natural n = 1; n < LN_ITERATIONS; ++n) {
    result += power((x-1)/(x+1), 2*n-1)/(2*n-1);
  }

  return 2*result;
}

// General logarithm -> loga(x)
MLDEF real ml_log(real base, real x) {
  return ml_ln(x) / ml_ln(base); // logb(x) = ln(x)/ln(b)
}

// Exponential function -- e^x
MLDEF real ml_exp(real x) {
  real result = 0;

  for(natural n = 0; n < EXP_ITERATIONS; ++n) {
    result += power(x, n) / (real)fact_norecurse(n);
  }

  return result;
}

// Used for exponential a^b -- allows for fractional and negative powers.
MLDEF real pow_exp(real base, real exponent) {
  return ml_exp( ml_ln(base) * exponent ); // a^b = e^(ln(a) * b)
}

// Gives an approximate square root of a number.
MLDEF real ml_sqrt(real x) {
  return pow_exp(x, 0.5);
}

#endif // ML_CORE_H
