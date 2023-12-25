/* The following values can be changed for performance or accuracy purposes:
 *  - INTEGRAL_ITERATIONS
 *  - SIN_ITERATIONS
 *  - LN_ITERATIONS
 *  - EXP_ITERATIONS
 *
 *  Note: Setting some of these values too high can result in related functions returning NAN or INFINITY due to accuracy issues.
 *  Changing these values waives your right to complain about performance.
 *  To change these values, define them BEFORE including this header file.
 */

#ifndef MATHLIB_H
#define MATHLIB_H

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

// Combinations formula
MLDEF natural nCr(natural n, natural r) {
  return fact(n) / (fact(r) * fact(n-r));
}

// Permutations formula
MLDEF natural nPr(natural n, natural r) {
  return fact(n) / fact(n-r);
}

// P(X=r) where X~B(n,p)
MLDEF real binomialProbability(natural n, natural r, real p) {
  return nCr(n, r) * power(p, r) * power(1-p, n-r);
}

// P(X<=r) where X~B(n,p)
MLDEF real binomialCumulativeProbability(natural n, natural r, real p) {
  real sum = 0;

  for(natural k = 0; k <= r; ++k) {
    sum += binomialProbability(n, k, p);
  }

  return sum;
}

MLDEF real binomialMean(natural n, real p) { return (real)n * p; }
MLDEF real binomialMedian(natural n, real p) { return ml_floor( binomialMean(n, p) ); }
MLDEF real binomialMode(natural n, real p) { return ml_floor( binomialMean(n+1, p) ); }
MLDEF real binomialVariance(natural n, real p) { return (real)n * p * (1-p); }
MLDEF real binomialStandardDeviation(natural n, real p) { return ml_sqrt( binomialVariance(n, p) ); }

MLDEF real arcsin(real x) {
  real sum = 0;

  for(natural n = 0; n < ARCSIN_ITERATIONS; ++n) {
    sum += ((real)fact(2*n) / (power(2.0, 2*n)*power((real)fact(2*n), 2))) * (power(x, 2*n+1)/(real)(2*n+1));
  }

  return sum;
}

MLDEF real arccos(real x) {
  return pi/2 - arcsin(x);
}

MLDEF real arctan(real x) {
  real sum = 0;

  if(x <= -1) {    
    for(natural n = 0; n < ARCTAN_ITERATIONS; ++n) {
      sum += power(-1, n) * (1/((2*n+1) * power(x, 2*n+1)));
    }

    sum = -pi/2 - sum;
  } else if(x <= 1) {
    for(natural n = 0; n < ARCTAN_ITERATIONS; ++n) {
      sum += power(-1, n) * (power(x, 2*n + 1)/(2*n+1));
    }
  } else {
    for(natural n = 0; n < ARCTAN_ITERATIONS; ++n) {
      sum += power(-1, n) * (1/((2*n+1) * power(x, 2*n+1)));
    }

    sum = pi/2 - sum;
  }
  
  return sum;
}

struct _complex {
  real realPart;
  real imaginaryPart;
};

typedef struct _complex complex;

MLDEF real argument(complex z) {
  return arctan(z.imaginaryPart / z.realPart);
}

MLDEF real modulus(complex z) {
  return ml_sqrt(z.realPart*z.realPart + z.imaginaryPart * z.imaginaryPart);
}

MLDEF complex conjugate(complex z) {
  return (complex) {
    .realPart = z.realPart,
    .imaginaryPart = -z.imaginaryPart
  };
}

MLDEF complex complex_sqrt(real x) {
  if(x < 0) {
    return (complex) {
      .realPart = 0,
      .imaginaryPart = ml_sqrt( ml_abs(x) ) 
    };
  } else {
    return (complex) {
      .realPart = ml_sqrt(x),
      .imaginaryPart = 0
    };
  }
}

MLDEF complex complex_add(complex a, complex b) {
  return (complex) {
    .realPart = a.realPart + b.realPart,
    .imaginaryPart = a.imaginaryPart + b.imaginaryPart
  };
}

MLDEF complex complex_sub(complex a, complex b) {
  return (complex) {
    .realPart = a.realPart - b.realPart,
    .imaginaryPart = a.imaginaryPart - b.imaginaryPart
  };
}

MLDEF complex complex_mul(complex a, complex b) {
  return (complex) {
    .realPart = a.realPart * b.realPart - a.imaginaryPart * b.imaginaryPart,
    .imaginaryPart = a.realPart * b.imaginaryPart + b.realPart * a.imaginaryPart
  };
}

MLDEF complex complex_div(complex a, complex b) {
  return (complex) {
    .realPart = (a.realPart * b.realPart + a.imaginaryPart * b.imaginaryPart)/(b.realPart * b.realPart + b.imaginaryPart * b.imaginaryPart),
    .imaginaryPart = (b.realPart * a.imaginaryPart - a.realPart * b.imaginaryPart)/(b.realPart * b.realPart + b.imaginaryPart * b.imaginaryPart)
  };
}

#endif // MATHLIB_H
