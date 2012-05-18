/*
 * fft_jp_2d_test.c
 *
 * Main test program for my fft_2d code.
 * This uses MPI for parallelization.
 *
 * Author: Jeffrey Picard
 */
/* Standard Headers */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

/* My headers */
#include "fft_jp.h"

/* Constants */
#define N 1024
#define M N*N
#define L 2*M_PI

#define DEBUG 1

/*
 * No arguments will run the normal fft_2d.
 * More than one argument will run the fft_2d_ret.
 */
int main( int argc, char **argv )
{
  MPI_Init( &argc, &argv );
  complex double *in = NULL, *out = NULL, *grid = NULL;
  FILE * fp = NULL;
  int rank = 0;

  in = calloc( N, sizeof(complex double) );
  if( !in )
    EXIT_WITH_PERROR("malloc failed in main: ");
  arr_fill_cosine( in, N, L );

  if( argc == 1 )
    fft_2d( N, in );
  else
  {
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if( rank == 0 )
    {
      out = malloc( M*sizeof(complex double) );
      grid = calloc( M, sizeof(complex double) );
      if( !out || !grid )
        EXIT_WITH_PERROR("malloc failed in main: ")
      arr_fill_grid_mpi( grid, M, 0, 1, 1 );
      fp = fopen("fft_jp_2d_ret.dat", "w");
      if( !fp )
        EXIT_WITH_PERROR("failed to open file in main: ");
    }

    fft_2d_ret( N, in, out );

    if( rank == 0 )
      write_data_arr( fp, out, grid, M, W_REAL );
  }

  free( in );
  free( out );
  free( grid );

  MPI_Finalize();
  return 0;
}
