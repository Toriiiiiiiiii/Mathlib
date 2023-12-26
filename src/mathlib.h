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

#include "mathlib_core.h"
#include "mathlib_discrete.h"
#include "mathlib_pure.h"
#include "mathlib_complex.h"
#include "mathlib_vector.h"

#endif // MATHLIB_H
