#ifndef ML_COMPLEX_H
#define ML_COMPLEX_H

#include "mathlib_core.h"
#include "mathlib_pure.h"

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

#endif // ML_COMPLEX_H
