#include <stdio.h>
#include "mathlib.h"

int main() {
  printf("ln(e) = %Lf\n", ml_ln(e));
  printf("log2(64) = %Lf\n", ml_log(2, 64));
  printf("2^3 = e^(ln(2) * 3) = %Lf (%Lf)\n", pow_exp(2, 3), power(2, 3));
  printf("sqrt(25) = %Lf\n", ml_sqrt(25));

  return 0;
}
