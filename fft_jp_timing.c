/*
 * fft_jp_timing.c
 *
 * Timing test program for my fft_jp code.
 *
 * Author: Jeffrey Picard
 */
/* Standard Headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* My headers */
#include "fft_jp.h"

/* Constants */
#define N 524288
#define M N
#define L 2*M_PI
#define MAX_TESTS 100

int main( int argc, char **argv )
{
  vector *in, *out; /* *grid;*/
  int num_tests = 0;
  int n = 2;
  /* Just in case, might be nice to have these */
  int points_arr[MAX_TESTS];
  double time_data[MAX_TESTS];
  double t_b, t_e;
  int i;
  if( argc == 2 )
  {
    num_tests = atoi( argv[1] );
    if( num_tests > MAX_TESTS )
      EXIT_WITH_ERROR("Error: max number of tests is %d\n", MAX_TESTS );
  }
  else
    EXIT_WITH_ERROR("Usage: fft_jp_test <num_pows_of_two>\n");
  FILE *fp = fopen("fft_jp_timing.dat", "w");
  if( !fp )
    EXIT_WITH_PERROR("file open failed in main: ")
  for( i = 0; i < num_tests; i++ )
  {
    t_b = WTime();
    create_vector( &in, n );
    create_vector( &out, n );
    /*create_vector( &grid, n );*/

    /*vec_fill_grid( grid, L );*/
    vec_fill_cosine( in, n, L );

    fft( in, out, n );

    /* write_data( fp, out, grid, N, W_REAL );*/

    destroy_vector( in );
    destroy_vector( out );
    t_e = WTime();
    points_arr[i] = n;
    time_data[i] = t_e - t_b;
    fprintf( fp, "%d %g\n", n, t_e - t_b );
    fprintf( stderr, "%d n = %d\n", i, n );
    /*free( grid );*/
    n *= 2;
  }

  fclose( fp );

  return 0;
}
