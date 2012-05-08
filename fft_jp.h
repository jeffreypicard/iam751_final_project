/*
 * fft_jp.h
 *
 * Header file for fft_jp.c
 *
 * Author: Jeffrey Picard (jpicardnh@gmail.com)
 */
#ifndef FFT_JP_H
#define FFT_JP_H

/* Headers */
#include <assert.h>
#include <complex.h>

/* Vector struct and macro */

typedef complex double val_type;

struct _vector {
  val_type *vals;
  int n;
} typedef vector;

#ifdef BOUNDSCHECK
#define VEC(v, i) (*({              \
  assert((i) >= 0 && (i) < (v)->n); \
  &((v)->vals[(i)]);                \
      })) 
#else
#define VEC(v, i) ((v)->vals[i])
#endif

#define EXIT_WITH_ERROR(...) {    \
  fprintf( stderr, __VA_ARGS__ ); \
  exit(-1);                       \
}

#define EXIT_WITH_PERROR(...) {   \
  fprintf( stderr, __VA_ARGS__ ); \
  perror( NULL );                 \
  exit(-1);                       \
}

/* Function Prototypes */
void fft( vector *, vector *, int, complex double );
void create_vector( vector **, int );
void destroy_vector( vector * );
void fill_grid_vec( vector *, val_type );

#endif
