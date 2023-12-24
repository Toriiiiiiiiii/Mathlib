/* Default number of iterations in an integral is 1,000,000.
 * To change this, define INTEGRAL_ITERATIONS as the desired number of iterations before including this file.
 */

#ifndef MATHLIB_H
#define MATHLIB_H

#ifndef INTEGRAL_ITERATIONS
  #define INTEGRAL_ITERATIONS 1000000
#endif // INTEGRAL_ITERATIONS

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

typedef double real;
typedef int integer;
typedef unsigned int natural;

typedef real (*function)(real n);

// Function that simply returns the input.
MLDEF real constant(real x) {return x;}

// Performs a discrete summation on a function -> Equivalent to Capital Sigma.
MLDEF real discrete_sum( function f, real start, real end, real step ) {
  real sum = 0;

  for(real n = start; n != end; n += step) {
    sum += CALL(f, n);
  }

  return sum;
}

// Performs a discrete multiplication on a function -> Equivalent to Capital Pi.
MLDEF real discrete_product( function f, real start, real end, real step ) {
  real sum = 1;

  for(real n = start; n != end; n += step) {
    sum *= CALL(f, n);
  }

  return sum;
}

// Uses recursion to calculate the factorial of a natural number.
MLDEF natural fact(natural n) {
  if(n == 0) return 1;

  return n * fact(n-1);
}

// Uses discrete_product to calculate the factorial of a natural number.
MLDEF natural fact_norecurse(natural n) {
  return discrete_product(constant, 1, n, 1);
}

// Gives the limit of a function as it approaches a value. Returns NAN if the limit is undefined.
MLDEF real lim(function f, real approach) {
  const real maxdiff = 0.0001f;

  real fh = CALL(f, approach + SMALL);
  real fh2 = CALL(f, approach - SMALL);

  if(fh - fh2 >= maxdiff) return NAN;

  return fh;
}

// Approximates a derivative at point x.
MLDEF real derivative(function f, real x) {
  real fx = CALL(f, x);
  real fh = CALL(f, x + SMALL);

  return (fh - fx)/(SMALL);
}

// Approximates a definite integral from a to b using the trapezium rule. 
MLDEF real integral(function f, real a, real b) {
  real area = 0;
  real x = a;

  real n = INTEGRAL_ITERATIONS;
  real step = (b - a)/n;

  while(x < b) {
    real x1 = x + step;
    area += CALL(f, x) + CALL(f, x1);

    x = x1;
  }

  return area * (step/2);
}

#endif // MATHLIB_H
