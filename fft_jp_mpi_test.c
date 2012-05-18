/*
 * fft_jp_mpi_test.c
 *
 * Main test program for my fft_jp code using MPI.
 *
 * FIXME
 *
 * I'm leaving this code the way it is becuase I got
 * everything working  here, but it is extremely messy.
 * The REAL main test for the parallel version of the
 * fft is in fft_jp_2d_test.c.
 *
 * FIXME
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
#define M N*N
#define L 2*M_PI

#define DEBUG 0
#define TAG   0x2357

int main( int argc, char **argv )
{
  MPI_Init(&argc, &argv);
  complex double **in = NULL, **out = NULL, *grid = NULL;
  int rank = 0, size = 0;
  int i = 0, j = 0;
  int chunk = N;
  int num_columns = 0;
  char fname[100];
  FILE *fp = NULL;
  complex double **send_mat = NULL;
  complex double **recv_mat = NULL;

  /* Make sure number of points is a power of 2 */
  assert( N % 2 == 0 );

  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );

  /* Make sure number of processes is a power of 2 */
  assert( size % 2 == 0 );

  num_columns = N / size;

#if DEBUG
  fprintf( stderr, "num_columns: %d\n", num_columns );
#endif

  /*
  in = malloc( num_columns*sizeof(complex double*) );
  out = malloc( num_columns*sizeof(complex double*) );
  grid = calloc( num_columns*chunk, sizeof(complex double) );
  if( !in || !out || !grid )
    EXIT_WITH_PERROR("malloc failed in main: ");

  for( i = 0; i < num_columns; i++ )
  {
    in[i] = calloc( chunk, sizeof(complex double) );
    out[i] = calloc( chunk, sizeof(complex double) );
    if( !in[i] || !out[i] )
      EXIT_WITH_PERROR("malloc failed in main: ");
    arr_fill_cosine( in[i], chunk, L );
  }*/

  in = malloc( num_columns*sizeof(complex double*) );
  out = malloc( num_columns*sizeof(complex double*) );
  grid = calloc( num_columns*chunk, sizeof(complex double) );
  if( !in || !out || !grid )
    EXIT_WITH_PERROR("malloc failed in main: ");

  in[0] = calloc( num_columns*chunk, sizeof(complex double) );
  out[0] = calloc( num_columns*chunk, sizeof(complex double) );

  if( !in[0] || !out[0] )
    EXIT_WITH_PERROR("malloc failed in main: ");

  for( i = 1; i < num_columns; i++ )
  {
    in[i] = in[i-1]+chunk;
    out[i] = out[i-1]+chunk;
  }

  for( i = 0; i < num_columns; i++ )
    arr_fill_cosine( in[i], chunk, L );
  arr_fill_grid_mpi( grid, chunk, rank, size, num_columns );

  /* Setup the matrices needed for alltoall communication */
  send_mat = malloc( size*sizeof(complex double*));
  //send_mat[0] = malloc( size*chunk*sizeof(complex double));
  send_mat[0] = malloc( size*chunk*num_columns*sizeof(complex double));
  for( i = 1; i < size; i++ )
    send_mat[i] = send_mat[i-1]+(chunk*num_columns);

  recv_mat = malloc( size*sizeof(complex double*));
  recv_mat[0] = malloc( size*chunk*num_columns*sizeof(complex double));
  for( i = 1; i < size; i++ )
    recv_mat[i] = recv_mat[i-1]+(chunk*num_columns);

#if DEBUG
  fprintf( stderr, "initialized matrices\n");
#endif

  /* phase 1, sqrt(N) fft */
  for( i = 0; i < num_columns; i++ )
    fft_mpi( in[i], out[i], chunk );

#if DEBUG
  fprintf( stderr, "finished first round of 1d fft's\n");
#endif

  /* Fills the matrix so all other processes get this same info */
  int k;
  for( i = 0; i < size; i++ )
    for( j = 0; j < num_columns; j++)
      for( k = 0; k < chunk; k++ )
        send_mat[i][(j*chunk)+k] = out[j][k];

#if DEBUG
  fprintf( stderr, "filled send_mat with data\n");
#endif

  /* phase 2, transpose */
  MPI_Alltoall( send_mat[0], chunk*num_columns, MPI_DOUBLE_COMPLEX, 
                recv_mat[0], chunk*num_columns, MPI_DOUBLE_COMPLEX,
                MPI_COMM_WORLD );

#if DEBUG
  fprintf( stderr, "finished MPI_Alltoall\n");
#endif

  /* 
   * FIXME There has to be a cleaner way to do this.
   * Until then here's an explanation.
   * We need to fill each column our process is responsible for. Hence
   * the outer loop.
   * We need to get infomation from every other process. Hence the middle
   * loop.
   * Each other process has num_columns and therefore we need num_columns 
   * elements from it. Hence the inner loop.
   * num_columns = chunk / size.
   * Therefore we may index through an entire column with 
   * n contained in [0...chunk] = [0...size*num_columns].
   * Hence the (i*num_columns)+k
   * Since each column in the recv_mat matrix has num_columns elements of
   * data we want, and this data is seperated by chunk number of elements 
   * and since we want the row in each other column which is equal to our 
   * column in the overall matrix, we want an index of (rank*num_columns)+j
   * into each column, therefore we index into the recv_mat matrix with 
   * (k*chunk)+(rank*num_columns)+j
   */
  for( j = 0; j < num_columns; j++ )
    for( i = 0; i < size; i++ )
      for( k = 0; k < num_columns; k++ )
        in[j][(i*num_columns)+k] = recv_mat[i][(k*chunk)+(rank*num_columns)+j];

#if DEBUG
  fprintf( stderr, "reloaded in data from recv_mat\n");
#endif

  /* phase 3, sqrt(N) fft */
  for( i = 0; i < num_columns; i++ )
    fft_mpi( in[i], out[i], chunk );

#if DEBUG
  fprintf( stderr, "finished second round of 1d fft's\n");
#endif


  /* phase 4, another transpose */
  /* Fills the matrix so all other processes get this same info */
  for( i = 0; i < size; i++ )
    for( j = 0; j < num_columns; j++)
      for( k = 0; k < chunk; k++ )
        send_mat[i][(j*chunk)+k] = out[j][k];
#if DEBUG
  fprintf( stderr, "filled send_mat with data second time\n");
#endif

  /* phase 5, retranspose */
  MPI_Alltoall( send_mat[0], chunk*num_columns, MPI_DOUBLE_COMPLEX, 
                recv_mat[0], chunk*num_columns, MPI_DOUBLE_COMPLEX,
                MPI_COMM_WORLD );

#if DEBUG
  fprintf( stderr, "finished MPI_Alltoall second time\n");
#endif

  /* FIXME same comment as above about this being hackish*/
  for( j = 0; j < num_columns; j++ )
    for( i = 0; i < size; i++ )
      for( k = 0; k < num_columns; k++ )
        in[j][(i*num_columns)+k] = recv_mat[i][(k*chunk)+(rank*num_columns)+j];

#if DEBUG
  fprintf( stderr, "reloaded in data from recv_mat second time\n");
#endif

  /*
  for( i = 0; i < num_columns; i++ )
  {
    sprintf(fname, "fft_jp_mpi%d%05d.dat", rank, i );
    fp = fopen(fname, "w");
    if( !fp )
      EXIT_WITH_PERROR("file open failed in main: ");

    write_data_arr( fp, in[i], (grid+(i*chunk)), chunk, W_REAL );
    fclose( fp );
  }*/
  sprintf(fname, "fft_jp_mpi%d.dat", rank );
  fp = fopen(fname, "w");
  if( !fp )
    EXIT_WITH_PERROR("file open failed in main: ");

  for( i = 0; i < num_columns; i++ )
  {
    write_data_arr( fp, in[i], (grid+(i*chunk)), chunk, W_REAL );
  }
  fclose( fp );


#if DEBUG
  fprintf( stderr, "finished writing to files.\n");
#endif
/*
  if( rank == 0 )
  {
    complex double *result = NULL;
    complex double *r_grid = NULL;
    result = malloc( size*num_columns*chunk*sizeof(complex double));
    r_grid = malloc( size*num_columns*chunk*sizeof(complex double) );
    if( !result || !r_grid )
      EXIT_WITH_PERROR("malloc failed in main: ");
    arr_fill_grid_mpi( r_grid, chunk, L, 0, 1, num_columns*size );
    k = 0;
    for( i = 0; i < num_columns; i++ )
      for( j = 0; j < chunk; j++ )
        result[k++] = in[i][j];
    for( i = 1; i < size; i++ )
      for( j = 0; j < num_columns; j++ )
      {
        //MPI_Request request;
        MPI_Recv( (result+(i*num_columns*chunk)+(j*chunk)), chunk, 
                   MPI_DOUBLE_COMPLEX, i, TAG, MPI_COMM_WORLD, 
                   MPI_STATUS_IGNORE );
      }
    
    fp = fopen("mpi_big.dat", "w");
    if( !fp )
      EXIT_WITH_PERROR("file open failed in main: ");
    write_data_arr( fp, result, r_grid, size*num_columns*chunk, W_REAL );
    fclose( fp );
  }
  else
  {
    for( i = 0; i < num_columns; i++ )
      MPI_Send( in[i], chunk, MPI_DOUBLE_COMPLEX, 0, TAG, MPI_COMM_WORLD );
  }*/

  free( in );
  free( out );
  free( grid );
  free( send_mat );
  free( recv_mat );

  MPI_Finalize();
  return 0;
}
