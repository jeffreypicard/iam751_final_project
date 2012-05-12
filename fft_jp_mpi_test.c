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
#define N 16
#define L 2*M_PI
#define W_REAL 0

int main( int argc, char **argv )
{
  MPI_Init(&argc, &argv);
  vector *in, *out, *grid;

  create_vector( &in, N );
  create_vector( &out, N );
  create_vector( &grid, N );

  vec_fill_grid( grid, L );
  vec_fill_cosine( in, L );

  int rank, size;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );

  char fname[100];
  sprintf(fname, "fft_jp_mpi%d", rank );
  FILE *fp = fopen(fname, "w");
  if( !fp )
    EXIT_WITH_PERROR("file open failed in main: ")

  int chunk = N/size;

  fprintf( stderr, "%d chunk: %d\n", rank, chunk );

  //int i_b = (N/size)*rank;
  //int i_e = (N/size)*(rank+1);

  vector *sub_in, *sub_out;
  create_vector( &sub_in, chunk );
  create_vector( &sub_out, chunk );
  int i, j;
  for( i = rank, j = 0; i < N; i += size, j++ )
  {
    VEC( sub_in, j ) = VEC( in, i );
  }

  write_vector( stderr, sub_in );
  write_vector( stderr, sub_out );
  /* phase 1, sqrt(N) fft */
  fft_mpi( sub_in, sub_out, chunk );

  /* phase 2, transpose */
  vector *sub_in2;
  vector *sub_out2;
  create_vector( &sub_in2, N );
  create_vector( &sub_out2, N );
  MPI_Alltoall( sub_out->vals, chunk, MPI_DOUBLE_COMPLEX, 
                sub_in2->vals, chunk, MPI_DOUBLE_COMPLEX,
                MPI_COMM_WORLD );

  /* phase 3, sqrt(N) fft */
  fft_mpi( sub_in2, sub_out2, chunk );

  write_data( fp, sub_out, grid, N/size, W_REAL );

  free( in );
  free( out );
  free( grid );

  fclose( fp );

  MPI_Finalize();
  return 0;
}
