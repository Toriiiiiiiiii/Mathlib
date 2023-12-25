/* Default number of iterations in an integral is 1,000,000.
 * To change this, define INTEGRAL_ITERATIONS as the desired number of iterations before including this file.
 */

#ifndef MATHLIB_H
#define MATHLIB_H

#ifndef INTEGRAL_ITERATIONS
  #define INTEGRAL_ITERATIONS 1000000
#endif // INTEGRAL_ITERATIONS

#ifndef SIN_ITERATIONS
  #define SIN_ITERATIONS 4
#endif // TAYLOR_ITERATIONS

#ifndef LN_ITERATIONS
  #define LN_ITERATIONS 1000
#endif

#ifndef EXP_ITERATIONS
  #define EXP_ITERATIONS 50
#endif

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

MLDEF real power(real x, natural n) {
  real result = 1;
  for(natural i = 0; i < n; ++i) {
    result *= x;
  }

  return result;
}

const real pi = 22.0 / 7.0;
const real e = 2.7182818284590452354;

// Function that simply returns the input.
MLDEF real constant(real x) {return x;}

MLDEF real ml_abs(real x) {return (x >= 0)? x : (x*-1);}

MLDEF real ml_floor(real x) { return (real)((natural)x); }
MLDEF real ml_ceil(real x) { 
  real fractionalPart = x - ml_floor(x);

  return (fractionalPart == 0)? x : (ml_floor(x)+1);
}

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
  natural result = 1;

  for(int x = 1; x < n+1; ++x) {
    result *= x;
  }

  return result;
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

// Converts degrees to radians
MLDEF real degrees_to_radians(real theta) {
  return theta * (pi / 180.0);
}

// Converts degrees to radians
MLDEF real radians_to_degrees(real theta) {
  return theta / (pi / 180.0);
}

// Sine trigonometric function.
MLDEF real ml_sin(real theta) {
  if(theta < 0) return -ml_sin( ml_abs(theta) );

  // For optimisation and accuracy purposes, we should only calculate values of theta for the range 0 - pi/2.
  // This allows for the number of iterations to be significantly reduced, and also improves the speed of the calculations.
  if(theta > pi/2) {
    natural period = ml_ceil( theta / pi );

    if(period % 2 == 0) {
      return -ml_sin( pi*period - theta );
    } else {
      return ml_sin( pi * period - theta );
    }
  }

  if(theta < 0.2 && theta > -0.2) return theta;

  real result = 0;
  for(natural n = 0; n < SIN_ITERATIONS; ++n) {
    natural n2 = (2*n) + 1;

    result += power(-1.0, n) * (power(theta, n2)/(real)fact_norecurse(n2));
  }

  return result;
}

// Cosine trigonometric function -> equivalent to -sin(x - pi/2)
MLDEF real ml_cos(real theta) {
  return -ml_sin( theta - pi/2 );
}

// Tanjent trigonometric function -> equivalent to sin(x)/cos(x)
MLDEF real ml_tan(real theta) {
  return ml_sin(theta) / ml_cos(theta);
}

// Natural logarithm -- logarithm with base e
MLDEF real ml_ln(real x) {
  real result = 0;

  for(natural n = 1; n < LN_ITERATIONS; ++n) {
    result += power((x-1)/(x+1), 2*n-1)/(2*n-1);
  }

  return 2*result;
}

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

#endif // MATHLIB_H
