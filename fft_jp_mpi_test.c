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
  //vector *in, *out, *grid;
  complex double *in, *out, *grid;
  complex double **in2, **out2, *grid2;
  //vector *sub_in, *sub_out;
  complex double *sub_in, *sub_out;
  int rank, size;
  int i, j;
  int chunk = N;
  int num_columns;
  char fname[100];
  FILE *fp = NULL;
  complex double **send_mat = NULL;
  complex double **recv_mat = NULL;
  //vector *sub_in2, *sub_out2;
  complex double *sub_in2, *sub_out2;
  //vector *result;
  complex double *result;

  /* Make sure number of points is a power of 2 */
  assert( N % 2 == 0 );

  in = calloc( N, sizeof(complex double) );
  out = calloc( N, sizeof(complex double) );
  grid = calloc( N, sizeof(complex double) );
  /*create_vector( &in, N );
  create_vector( &out, N );
  create_vector( &grid, M );*/

  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );

  /* Make sure number of processes is a power of 2 */
  assert( size % 2 == 0 );

  num_columns = M / size;

  in2 = malloc( num_columns*sizeof(complex double*) );
  out2 = malloc( num_columns*sizeof(complex double*) );
  grid2 = calloc( num_columns*chunk, sizeof(complex double) );
  if( !in2 || !out2 || !grid2 )
    EXIT_WITH_PERROR("malloc failed in main: ");

  for( i = 0; i < num_columns; i++ )
  {
    in2[i] = calloc( chunk, sizeof(complex double) );
    out2[i] = calloc( chunk, sizeof(complex double) );
    if( !in2[i] || !out2[i] )
      EXIT_WITH_PERROR("malloc failed in main: ");
    arr_fill_cosine( in2[i], chunk, L );
  }
  arr_fill_grid_mpi( grid2, chunk*num_columns, L, rank, size );

  arr_fill_grid_mpi( grid, N, L, rank, size );
  arr_fill_cosine( in, N, L );



  fprintf( stderr, "%d chunk: %d\n", rank, chunk );

  //int i_b = (N/size)*rank;
  //int i_e = (N/size)*(rank+1);

  /*
  sub_in = calloc( chunk, sizeof(complex double) );
  sub_out = calloc( chunk, sizeof(complex double) );
  create_vector( &sub_in, chunk );
  create_vector( &sub_out, chunk );
  for( i = rank, j = 0; i < N; i += size, j++ )
  {
    VEC( sub_in, j ) = VEC( in, i );
  }
  for( i = 0; i < N; i++ )
    sub_in[i] = in[i];*/
    //VEC( sub_in, i ) = VEC( in, i );

  /* Setup the matrices needed for alltoall communication */
  send_mat = malloc( size*sizeof(complex double*));
  //send_mat[0] = malloc( size*chunk*sizeof(complex double));
  send_mat[0] = malloc( size*chunk*num_columns*sizeof(complex double));
  for( i = 1; i < size; i++ )
    send_mat[i] = send_mat[i-1]+(chunk*num_columns);

  recv_mat = malloc( size*sizeof(complex double*));
  //recv_mat[0] = malloc( size*chunk*sizeof(complex double));
  recv_mat[0] = malloc( size*chunk*num_columns*sizeof(complex double));
  for( i = 1; i < size; i++ )
    recv_mat[i] = recv_mat[i-1]+(chunk*num_columns);

  //write_arr( stderr, sub_in, N );
  //write_vector( stderr, sub_out );
  /* phase 1, sqrt(N) fft */
  //fft_mpi( sub_in, sub_out, chunk );
  for( i = 0; i < num_columns; i++ )
    fft_mpi( in2[i], out2[i], chunk );

  /* Fills all the matrices so all other processes get this same info */
  //int e_j = chunk*num_columns;
  int k;
  for( i = 0; i < size; i++ )
    for( j = 0; j < num_columns; j++)
      for( k = 0; k < chunk; k++ )
        send_mat[i][(j*chunk)+k] = out2[j][k];

  /* phase 2, transpose */
  //sub_in2 = calloc( chunk, sizeof(complex double) );
  //sub_out2 = calloc( chunk, sizeof(complex double) );
  /*create_vector( &sub_in2, chunk );
  create_vector( &sub_out2, chunk );*/
  MPI_Alltoall( send_mat[0], chunk, MPI_DOUBLE_COMPLEX, 
                recv_mat[0], chunk, MPI_DOUBLE_COMPLEX,
                MPI_COMM_WORLD );
  /* We want our rank's column for each row for a matrix tranpose. 
   * Again this assumes there are sqrt(n) proccess running */
  /*
  for( i = 0; i < chunk; i++ )
    sub_in2[i] = recv_mat[i][rank];
  write_arr( stderr, sub_in2, N );*/
  for( i = 0; i < size; i++ )
    for( j = 0; j < num_columns; j++ );
      for( k = 0; k < chunk; k++ )
        in2[(i*size)+j][k] = recv_mat[i][(j*chunk)+k];
  /* phase 3, sqrt(N) fft */
  for( i = 0; i < num_columns; i++ )
    fft_mpi( in2[i], out2[i], chunk );

  /* phase 4, another transpose */
  /* Fills all the matrices so all other processes get this same info */
  for( i = 0; i < size; i++ )
    for( j = 0; j < num_columns; j++)
      for( k = 0; k < chunk; k++ )
        send_mat[i][(j*chunk)+k] = out2[j][k];
  /*
  for( i = 0; i < size; i++ )
    for( j = 0; j < chunk; j++)
      send_mat[i][j] = sub_out2[j];*/

  //result = calloc( chunk, sizeof(complex double) );
  /*create_vector( &result, chunk );*/
  MPI_Alltoall( send_mat[0], chunk, MPI_DOUBLE_COMPLEX, 
                recv_mat[0], chunk, MPI_DOUBLE_COMPLEX,
                MPI_COMM_WORLD );
  /*
  for( i = 0; i < chunk; i++ )
    result[i] = recv_mat[i][rank];*/

  for( i = 0; i < size; i++ )
    for( j = 0; j < num_columns; j++ );
      for( k = 0; k < chunk; k++ )
        in2[(i*size)+j][k] = recv_mat[i][(j*chunk)+k];

  for( i = 0; i < num_columns; i++ )
  {
    sprintf(fname, "fft_jp_mpi%d%d.dat", rank, i );
    fp = fopen(fname, "w");
    if( !fp )
      EXIT_WITH_PERROR("file open failed in main: ")

    write_data_arr( fp, in2[i], grid2, chunk, W_REAL );
  }

  /*
  free( in );
  free( out );
  free( sub_in );
  free( sub_out );
  free( sub_in2 );
  free( sub_out2 );
  free( result );
  free( grid );*/
  free( in2 );
  free( out2 );
  free( send_mat );
  free( recv_mat );

  fclose( fp );

  MPI_Finalize();
  return 0;
}
