#include <stdio.h>
#include "mathlib.h"

real f(real x) {
  return x*x*x;
}

real reciprocal(real x) {
  return 0 / x;
}

int main() {
  printf("10! = %lu (%lu)\n", fact(10), fact_norecurse(10));

  for(natural theta = 0; theta < 360; theta += 45) {
    printf("sin(%lu) = %Lf\n", theta, ml_sin( degrees_to_radians((real)theta) ));
  }

  return 0;
}
