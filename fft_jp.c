/*
 * fft_jp.c
 *
 * A parallel fft implementation written as a final project
 * for IAM751 Spring 2012.
 *
 * Author: Jeffrey Picard (jpicardnh@gmail.com)
 */

/* Standard headers */
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <mpi.h>

/* My headers */
#include "fft_jp.h"

/*
 * fft
 *
 * Takes an input vector, output vector and a number of points
 * Performs a Fast Fourier Transform on the input vector
 * and stores it in the output vector.
 */
void fft( vector *in, vector *out, int n )
{
  int k;
  if( n == 1 )
    VEC( out, 0 ) = VEC( in, 0 );
  else
  {
    complex double wn = cexp( -2*M_PI*I / n );
    vector *x, *y, *p, *q;
    create_vector( &x, n/2 );
    create_vector( &y, n/2 );
    create_vector( &p, n/2 );
    create_vector( &q, n/2 );
    for( k = 0; k < n/2; k++ )
    {
      VEC( x, k ) = VEC( in, 2*k );
      VEC( y, k ) = VEC( in, 2*k + 1 );
    }
    fft( x, p, n/2 );
    fft( y, q, n/2 );
    for( k = 0; k < n; k++ )
    {
      VEC( out, k ) = VEC( p, k % (n/2) ) + cpow(wn,k)*VEC( q, k % (n/2) );
    }
    destroy_vector( x );
    destroy_vector( y );
    destroy_vector( p );
    destroy_vector( q );
  }
}

/*
 * fft_2d
 *
 * This function takes a number of points and performs
 * a parallel 2d fft on the input, which is assumed to
 * be in row major order. This function returns no results
 * but write the output to a file.
 */
void fft_2d( const int n, complex double *data )
{
  complex double **in = NULL, **out = NULL, *grid = NULL;
  int rank = 0, size = 0;
  int i = 0, j = 0;
  int chunk = n;
  int num_columns = 0;
  char fname[100];
  FILE *fp = NULL;
  complex double **send_mat = NULL;
  complex double **recv_mat = NULL;

  /* Make sure number of points is a power of 2 */
  assert( n % 2 == 0 );

  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );

  /* Make sure number of processes is divisible 2 */
  assert( size % 2 == 0 );

  num_columns = n / size;

#if DEBUG
  fprintf( stderr, "num_columns: %d\n", num_columns );
#endif

  in = malloc( num_columns*sizeof(complex double*) );
  out = malloc( num_columns*sizeof(complex double*) );
  grid = calloc( num_columns*chunk, sizeof(complex double) );
  if( !in || !out || !grid )
    EXIT_WITH_PERROR("malloc failed in fft_2d: ");

  in[0] = calloc( num_columns*chunk, sizeof(complex double) );
  out[0] = calloc( num_columns*chunk, sizeof(complex double) );

  if( !in[0] || !out[0] )
    EXIT_WITH_PERROR("malloc failed in fft_2d: ");

  for( i = 1; i < num_columns; i++ )
  {
    in[i] = in[i-1]+chunk;
    out[i] = out[i-1]+chunk;
  }

  /* copy over input data */
  for( i = 0; i < num_columns; i++ )
    for( j = 0; j < chunk; j++ )
      in[i][j] = data[j];

  arr_fill_grid_mpi( grid, chunk, rank, size, num_columns );

  /* Setup the matrices needed for alltoall communication */
  send_mat = malloc( size*sizeof(complex double*));
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

  /* Write data to files */
  sprintf(fname, "fft_jp_mpi%d.dat", rank );
  fp = fopen(fname, "w");
  if( !fp )
    EXIT_WITH_PERROR("file open failed in fft_2d: ");

  for( i = 0; i < num_columns; i++ )
  {
    write_data_arr( fp, in[i], (grid+(i*chunk)), chunk, W_REAL );
  }
  fclose( fp );

  /* Cleanup memory */
  free( in[0] );
  free( in );
  free( out[0] );
  free( out );
  free( recv_mat[0] );
  free( recv_mat );
  free( send_mat[0] );
  free( send_mat );
  free( grid );
}

/*
 * fft_2d_ret
 *
 * This function takes a number of points and performs
 * a parallel 2d fft on the input, which is assumed to
 * be in row major order. This function returns the result
 * in the out vector on the process with rank 0.
 *
 * WARNING: This is much slower due to the extra MPI
 * communication needed at the end!
 */
void fft_2d_ret( const int n, complex double *data, complex double *result )
{
  complex double **in = NULL, **out = NULL, *grid = NULL;
  int rank = 0, size = 0;
  int i = 0, j = 0;
  int chunk = n;
  int num_columns = 0;
  complex double **send_mat = NULL;
  complex double **recv_mat = NULL;
  const int TAG = 0x2357;

  /* Make sure number of points is a power of 2 */
  assert( n % 2 == 0 );

  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );

  /* Make sure number of processes is divisible of 2 */
  assert( size % 2 == 0 );

  num_columns = n / size;

#if DEBUG
  fprintf( stderr, "num_columns: %d\n", num_columns );
#endif

  in = malloc( num_columns*sizeof(complex double*) );
  out = malloc( num_columns*sizeof(complex double*) );
  grid = calloc( num_columns*chunk, sizeof(complex double) );
  if( !in || !out || !grid )
    EXIT_WITH_PERROR("malloc failed in fft_2d_ret: ");

  in[0] = calloc( num_columns*chunk, sizeof(complex double) );
  out[0] = calloc( num_columns*chunk, sizeof(complex double) );

  if( !in[0] || !out[0] )
    EXIT_WITH_PERROR("malloc failed in fft_2d_ret: ");

  for( i = 1; i < num_columns; i++ )
  {
    in[i] = in[i-1]+chunk;
    out[i] = out[i-1]+chunk;
  }

  /* copy over input data */
  for( i = 0; i < num_columns; i++ )
    for( j = 0; j < chunk; j++ )
      in[i][j] = data[j];

  arr_fill_grid_mpi( grid, chunk, rank, size, num_columns );

  /* Setup the matrices needed for alltoall communication */
  send_mat = malloc( size*sizeof(complex double*));
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

  /* Communicate the total result to the result array given to proccess 0 */
  if( rank == 0 )
  {
    k = 0;
    for( i = 0; i < num_columns; i++ )
      for( j = 0; j < chunk; j++ )
        result[k++] = in[i][j];

    for( i = 1; i < size; i++ )
      for( j = 0; j < num_columns; j++ )
      {
        MPI_Recv( (result+(i*num_columns*chunk)+(j*chunk)), chunk, 
                   MPI_DOUBLE_COMPLEX, i, TAG, MPI_COMM_WORLD, 
                   MPI_STATUS_IGNORE );
      }
  }
  else
  {
    for( i = 0; i < num_columns; i++ )
      MPI_Send( in[i], chunk, MPI_DOUBLE_COMPLEX, 0, TAG, MPI_COMM_WORLD );
  }

  /* Cleanup memory */
  free( in[0] );
  free( in );
  free( out[0] );
  free( out );
  free( recv_mat[0] );
  free( recv_mat );
  free( send_mat[0] );
  free( send_mat );
  free( grid );
}

/*
 * fft_mpi
 *
 * Takes an input vector, output vector and a number of points
 * Performs a Fast Fourier Transform on the input vector
 * and stores it in the output vector.
 *
 * Only difference between this and the fft function is that
 * this uses complex double arrays instead of the VEC struct.
 */
void fft_mpi( complex double *in, complex double *out, int n )
{
  int k;
  if( n == 1 )
    out[0] = in[0];
  else
  {
    complex double wn = cexp( -2*M_PI*I / n );
    complex double *x, *y, *p, *q;
    x = calloc( n/2, sizeof(complex double) );
    y = calloc( n/2, sizeof(complex double) );
    p = calloc( n/2, sizeof(complex double) );
    q = calloc( n/2, sizeof(complex double) );
    if( !x || !y || !p || !q )
      EXIT_WITH_PERROR("malloc failed in fft_mpi: ");
    for( k = 0; k < n/2; k++ )
    {
      x[k] = in[2*k];
      y[k] = in[2*k + 1];
    }
    fft_mpi( x, p, n/2 );
    fft_mpi( y, q, n/2 );
    for( k = 0; k < n; k++ )
      out[k] = p[ k % (n/2) ] + cpow(wn,k)*q[k % (n/2) ];

    free( x );
    free( y );
    free( p );
    free( q );
  }
}

/*
 * fft_bin
 *
 * Takes an input vector, output vector and a number of points
 * Performs a Fast Fourier Transform on the input vector
 * and stores it in the output vector.
 * 
 * Uses the binary exchange algorithm
 */
void fft_bin( vector *in, vector *out, int n )
{
  /*TODO*/
  /*Implement this at some point*/
}

/*
 * vec_fill_sine
 *
 * Takes a vector and a scalar and fills the vector
 * with sine values.
 */
void vec_fill_sine( vector *v, val_type scale )
{
  int i;
  for( i = 0; i < v->n; i++ )
    VEC( v, i ) = sin( i * scale / v->n );
}

/*
 * vec_fill_cosine
 *
 * Takes a vector and a scalar and fills the vector
 * with cosine values.
 */
void vec_fill_cosine( vector *v, const int n, val_type scale )
{
  int i;
  for( i = 0; i < n; i++ )
    VEC( v, i ) = cos( i * scale / n );
}

/*
 * arr_fill_cosine
 *
 * Takes an array, an int and a scalar and fills the vector
 * with cosine values.
 */
void arr_fill_cosine( complex double *v, const int n, complex double scale )
{
  int i;
  for( i = 0; i < n; i++ )
    v[i] = cos( i * scale / n );
}

/*
 * arr_fill_sine
 *
 * Takes an array, an int and a scalar and fills the vector
 * with sine values.
 */
void arr_fill_sine( complex double *v, const int n, complex double scale )
{
  int i;
  for( i = 0; i < n; i++ )
    v[i] = sin( i * scale / n );
}

/*
 * vec_fill_grid
 *
 * Takes a vector and a scalar and fills the vector
 * to be an evenly spaced grid from 0 ... scale.
 */
void vec_fill_grid( vector *v, val_type scale )
{
  int i;
  for( i = 0; i < v->n; i++ )
    VEC( v, i ) = (double)i / (double)v->n;
}

/*
 * vec_fill_grid_mpi
 *
 * Takes a vector and a scalar and fills the vector
 * to be an evenly spaced grid from 0 ... scale.
 * Scaling with the rank from MPI.
 */
void vec_fill_grid_mpi( vector *v, val_type scale, int rank, int size )
{
  int i;
  for( i = 0; i < v->n; i++ )
    VEC( v, i ) = (double)((rank*size)+i) / (double)v->n;
    //VEC( v, i ) = ((double)i / (double)v->n) * (double)(rank+1);
}

/*
 * arr_fill_grid_mpi
 *
 * Takes an array, the number of points per chunk, the rank
 * of the process calling, the total number of proccesses and
 * the number of columns the process is responsible for and
 * fills a grid array.
 */
void arr_fill_grid_mpi( complex double *v, const int n, const int rank, 
                        const int size, const int num_cols )
{
  int i;
  int e_i = num_cols*n;
  for( i = 0; i < e_i; i++ )
    v[i] = (double)((rank*e_i)+i) / (double)(n*num_cols*size);
}

/*
 * create_vector
 *
 * Takes a pointer to a pointer to a vector and a size.
 * Allocates space for the vector struct and the number
 * of specified values.
 */
void create_vector( vector **v, int n )
{
  (*v) = malloc( sizeof(vector) );
  if( !*v )
    EXIT_WITH_PERROR("malloc failed in create_vector: ")
  (*v)->vals = calloc( n, sizeof(val_type) );
  if( !(*v)->vals )
    EXIT_WITH_PERROR("malloc failed in create_vector: ")
  (*v)->n = n;
}

/*
 * destroy_vector
 *
 * Takes a vector pointer and frees all of its
 * associated memory, including the struct.
 */
void destroy_vector( vector *v )
{
  if( !v ) return;
  free( v->vals );
  free( v );
}

/*
 * write_vector
 *
 * Takes a FILE*, and a vector and writes it 
 */
void write_vector( FILE *fp, vector *v )
{
  int i;
  fprintf( fp, "\n" );
  for( i = 0; i < v->n; i++ )
    fprintf( fp, "%g ", (double)VEC( v, i ) );
  fprintf( fp, "\n" );
}

/*
 * write_arr
 *
 * Takes a FILE*, an array and an int and writes it 
 */
void write_arr( FILE *fp, complex double *v, const int n )
{
  int i;
  fprintf( fp, "\n" );
  for( i = 0; i < n; i++ )
    fprintf( fp, "%g ", (double)v[i] );
  fprintf( fp, "\n" );
}

/*
 * write_data
 *
 * Takes a FILE*, a data vector, a grid vector 
 * a size and a mode.
 * mode 0: write real values
 * mode 1: write complex values
 * mode 2: write both real and comples values
 */
void write_data( FILE *fp, vector *v, vector *grid, 
                      const int n, int mode )
{
  int i;
  if( mode == W_REAL )
    for( i = 0; i < n; i++ )
      fprintf( fp, "%g %g\n", creal(VEC( grid, i )), creal(VEC( v, i )) );
  else if( mode == W_COMP )
    for( i = 0; i < n; i++ )
      fprintf( fp, "%g %g\n", creal(VEC( grid, i )), cimag(VEC( v, i )) );
  else if( mode == W_BOTH )
    for( i = 0; i < n; i++ )
      fprintf( fp, "%g %g+i%g\n", creal(VEC( grid, i )), creal(VEC( v, i )),
                                                         cimag(VEC( v, i )));
  else
    EXIT_WITH_ERROR("Error: In write_data_arr, unkown mode.\n");
}

/*
 * write_data_arr
 *
 * Takes a FILE*, a data array, a grid array
 * a size and a mode.
 * mode 0: write real values
 * mode 1: write complex values
 * mode 2: write both real and comples values
 */
void write_data_arr( FILE *fp, complex double *v, complex double *grid, 
                      const int n, const int mode )
{
  int i;
  if( mode == W_REAL )
    for( i = 0; i < n; i++ )
      fprintf( fp, "%g %g\n", creal(grid[i]), creal(v[i]) );
  else if( mode == W_COMP )
    for( i = 0; i < n; i++ )
      fprintf( fp, "%g %g\n", creal(grid[i]), cimag(v[i]) );
  else if( mode == W_BOTH )
    for( i = 0; i < n; i++ )
      fprintf( fp, "%g %g+i%g\n", creal(grid[i]), creal(v[i]), cimag(v[i]) );
  else
    EXIT_WITH_ERROR("Error: In write_data_arr, unknown mode.\n");

}
