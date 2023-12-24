#include <stdio.h>
#include "mathlib.h"

real f(real x) {
  return x;
}

real reciprocal(real x) {
  return 0 / x;
}

int main() {
  real limit = lim(reciprocal, 0);

  if(ISINF(limit)) {
    printf("Limit of 0/x as x approaches 0 = infinity\n");
  } else if(ISNAN(limit)) {
    printf("The limit of 0/x as x approaches 0 is not defined.\n");
  } else {
    printf("Limit of 0/x as x approaches 0 = %lf\n", limit);
  }

  return 0;
}
