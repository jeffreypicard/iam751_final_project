/*
 * fft_jp_mpi_test.c
 *
 * Main test program for my fft_jp code using MPI.
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
#define N 4
#define M N*N
#define L 2*M_PI
#define W_REAL 0

/* 
 * FIXME 
 *
 * This is just trying to get MPI_Alltoall working with the fft
 * algorithm and get the parallel version actually getting a correct
 * answer. This currently assumes there are sqrt(N) processes running
 *
 * FIXME
 */
int main( int argc, char **argv )
{
  MPI_Init(&argc, &argv);
  vector *in, *out, *grid;

  create_vector( &in, N );
  create_vector( &out, N );
  create_vector( &grid, M );


  int rank, size;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );

  vec_fill_grid_mpi( grid, L, rank, size );
  vec_fill_cosine( in, N, L );


  char fname[100];
  sprintf(fname, "fft_jp_mpi%d.dat", rank );
  FILE *fp = fopen(fname, "w");
  if( !fp )
    EXIT_WITH_PERROR("file open failed in main: ")

  int chunk = N;

  fprintf( stderr, "%d chunk: %d\n", rank, chunk );

  //int i_b = (N/size)*rank;
  //int i_e = (N/size)*(rank+1);

  vector *sub_in, *sub_out;
  create_vector( &sub_in, chunk );
  create_vector( &sub_out, chunk );
  int i, j;
  /*
  for( i = rank, j = 0; i < N; i += size, j++ )
  {
    VEC( sub_in, j ) = VEC( in, i );
  }*/
  for( i = 0; i < N; i++ )
    VEC( sub_in, i ) = VEC( in, i );

  /* Setup the matrices needed for alltoall communication */
  complex double **send_mat = malloc( size*sizeof(complex double*));
  send_mat[0] = malloc( size*chunk*sizeof(complex double));
  for( i = 1; i < size; i++ )
    send_mat[i] = send_mat[i-1]+chunk;

  complex double **recv_mat = malloc( size*sizeof(complex double*));
  recv_mat[0] = malloc( size*chunk*sizeof(complex double));
  for( i = 1; i < size; i++ )
    recv_mat[i] = recv_mat[i-1]+chunk;

  write_vector( stderr, sub_in );
  //write_vector( stderr, sub_out );
  /* phase 1, sqrt(N) fft */
  fft_mpi( sub_in, sub_out, chunk );

  /* Fills all the matrices so all other processes get this same info */
  for( i = 0; i < size; i++ )
    for( j = 0; j < chunk; j++)
      send_mat[i][j] = sub_out->vals[j];

  /* phase 2, transpose */
  vector *sub_in2;
  vector *sub_out2;
  create_vector( &sub_in2, chunk );
  create_vector( &sub_out2, chunk );
  MPI_Alltoall( send_mat[0], chunk, MPI_DOUBLE_COMPLEX, 
                recv_mat[0], chunk, MPI_DOUBLE_COMPLEX,
                MPI_COMM_WORLD );
  /* We want our rank's column for each row for a matrix tranpose. 
   * Again this assumes there are sqrt(n) proccess running */
  for( i = 0; i < chunk; i++ )
    sub_in2->vals[i] = recv_mat[i][rank];
  write_vector( stderr, sub_in2 );
  /* phase 3, sqrt(N) fft */
  fft_mpi( sub_in2, sub_out2, chunk );

  /* phase 4, another transpose */
  /* Fills all the matrices so all other processes get this same info */
  for( i = 0; i < size; i++ )
    for( j = 0; j < chunk; j++)
      send_mat[i][j] = sub_out2->vals[j];
  vector *result;
  create_vector( &result, chunk );
  MPI_Alltoall( send_mat[0], chunk, MPI_DOUBLE_COMPLEX, 
                recv_mat[0], chunk, MPI_DOUBLE_COMPLEX,
                MPI_COMM_WORLD );
  for( i = 0; i < chunk; i++ )
    result->vals[i] = recv_mat[i][rank];

  write_data( fp, result, grid, N, W_REAL );

  free( in );
  free( out );
  free( grid );

  fclose( fp );

  MPI_Finalize();
  return 0;
}
