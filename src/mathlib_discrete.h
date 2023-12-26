#ifndef ML_DISCRETE_H
#define ML_DISCRETE_H

#include "mathlib_core.h"

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

#endif // ML_DISCRETE_H
