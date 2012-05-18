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

/* 
 * Originally these vector structs and macros
 * were used but I scrapped them when trying
 * to get the 2d fft working.
 */
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
void fft( vector *, vector *, int );
void fft_mpi( complex double *, complex double *, int );
void create_vector( vector **, int );
void destroy_vector( vector * );
void vec_fill_grid( vector *, val_type );
void vec_fill_grid_mpi( vector *, val_type, int, int );
void arr_fill_grid_mpi( complex double *, const int, 
                        const int, const int, const int );
void vec_fill_sine( vector *, val_type );
void arr_fill_sine( complex double *, const int, complex double );
void vec_fill_cosine( vector *, const int, val_type );
void arr_fill_cosine( complex double *, const int, complex double );

/* 
 * These constant go with the write_data functions and change 
 * their behavior when passed in as the last argument.
 */
#define W_REAL 0
#define W_COMP 1
#define W_BOTH 2
void write_data( FILE *, vector *, vector *, const int , int );
void write_data_arr( FILE *, complex double *, complex double *, 
                     const int, const int );
void write_vector( FILE *, vector * );
void write_arr( FILE *, complex double *, const int );

void fft_2d( const int, complex double *);
void fft_2d_ret( const int, complex double *, complex double *);

#include <sys/time.h>

static inline double
WTime( void )
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec / 1e6;
}

#endif
