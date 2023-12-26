#ifndef ML_PURE_H
#define ML_PURE_H

#include "mathlib_core.h"
#include "mathlib_discrete.h"

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

#endif // ML_PURE_H
