#ifndef ML_VECTOR_H
#define ML_VECTOR_H

#include "mathlib_core.h"

struct _vector {
  real dimensions;
  real *values;
};

typedef struct _vector vector;

#define VECTOR_UNDEFINED (vector){0, NULL}

// Creates a vector structure.
// Due to limitations, dimensions must be between 1 and 4 inclusive.
// The returned structure, if not VECTOR_UNDEFINED, should be destroyed using destroyVector() when it is no longer in use.
MLDEF vector createVector(real dimensions) {
  vector result = {0};

  if(dimensions < 1 || dimensions > 4) return VECTOR_UNDEFINED;

  result.dimensions = dimensions;
  result.values = (real *)malloc( dimensions * sizeof(real) );

  return result;
}

// Destroys a vector structure.
MLDEF void destroyVector(vector *v) {
  v->dimensions = 0;
  free(v->values);
}

MLDEF real getVectorComponent(vector v, natural index) {
  if(index >= v.dimensions) return NAN;

  return v.values[index];
}

MLDEF void setVectorComponent(vector *v, natural index, real value) {
  if(index >= v->dimensions) return;

  v->values[index] = value;
}

MLDEF vector vector_add(vector a, vector b) {
  if(a.dimensions != b.dimensions) return VECTOR_UNDEFINED;

  vector result = createVector(a.dimensions);
  for(natural index = 0; index < a.dimensions; ++index) {
    setVectorComponent(&result, index, 
        getVectorComponent(a, index) + getVectorComponent(b, index));
  }

  return result;
}

MLDEF vector vector_sub(vector a, vector b) {
  if(a.dimensions != b.dimensions) return VECTOR_UNDEFINED;

  vector result = createVector(a.dimensions);
  for(natural index = 0; index < a.dimensions; ++index) {
    setVectorComponent(&result, index, 
        getVectorComponent(a, index) - getVectorComponent(b, index));
  }

  return result;
}

MLDEF vector vector_scalar_multiply(vector a, real b) {
  vector result = createVector(a.dimensions);

  for(natural index = 0; index < a.dimensions; ++index) {
    setVectorComponent(&result, index, 
        getVectorComponent(a, index) * b);
  }

  return result;
}

MLDEF real vector_dot(vector a, vector b) {
  if(a.dimensions != b.dimensions) return NAN;

  real result = 0;
  for(natural index = 0; index < a.dimensions; ++index) {
    result += getVectorComponent(a, index) * getVectorComponent(b, index);
  }

  return result;
}

#endif // ML_VECTOR_H
