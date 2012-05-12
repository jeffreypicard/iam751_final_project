/*
 * fft_jp_mpi_test.c
 *
 * Main test program for my fft_jp code.
 *
 * Author: Jeffrey Picard
 */
/* Standard Headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

/* My headers */
#include "fft_jp.h"

/* Constants */
#define N 1024
#define L 2*M_PI
#define W_REAL 0

int main( int argc, char **argv )
{
  vector *in, *out, *grid;
  FILE *fp = fopen("fft_jp.dat", "w");
  if( !fp )
    EXIT_WITH_PERROR("file open failed in main: ")

  MPI_Init(&argc, &argv);

  create_vector( &in, N );
  create_vector( &out, N );
  create_vector( &grid, N );

  vec_fill_grid( grid, L );
  vec_fill_cosine( in, L );

  fft_mpi( in, out, N );

  write_data( fp, out, grid, N, W_REAL );

  free( in );
  free( out );
  free( grid );

  fclose( fp );

  MPI_Finalize();
  return 0;
}
