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
  printf("%Lf\n", ml_ceil(3.4f));

  for(natural theta = 0; theta <= 360; theta += 15) {
    printf("sin(%lu) = %Lf, cos(%lu) = %Lf, tan(%lu) = %Lf\n", theta, ml_sin( degrees_to_radians((real)theta) ), 
                                                               theta, ml_cos( degrees_to_radians((real)theta) ),
                                                               theta, ml_tan( degrees_to_radians((real)theta) ));
  }

  return 0;
}
