#include <stdio.h>
#include "mathlib.h"

int main() {
  printf("ln(e) = %Lf\n", ml_ln(e));
  printf("2^3 = e^(ln(2) * 3) = %Lf (%Lf)\n", pow_exp(2, 3), power(2, 3));

  return 0;
}
