/*------------------------------------------------------------
 File       : sparse2.c
 Created    : Thu Aug  7 16:10:33 2008
 Author     : Jostein R. Natvig <Jostein.R.Natvig@sintef.no>
 Revision   : $Id: sparse.c 186 2008-10-29 16:52:37Z jrn $


 Description


 Changes

------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <SuiteSparse/umfpack.h>

#include "lu.h"
#include "sparse.h"


#if MEXFILE
#include <mex.h>
#define print   mexPrintf
#define malloc  mxMalloc
#define calloc  mxCalloc
#define realloc mxRealloc
#define free    mxFree
#else

#endif


static const double EPS = 1e-14;
static void sort(const int n, int *row, int *col, double *value);

/*
 * Input arguments:
 *
 *     nrows    -
 *
 *     ncols    -
 *
 *     ndata    -
 *
 * Output:
 *
 *     return   - Sparse nrows x ncols matrix witgh room for ndata nonzeros.
 *
 */
/*---------------------------------------------------------------------------*/
sparse_t *sparse_alloc(int nrows, int ncols, int ndata)
/*---------------------------------------------------------------------------*/
{
  assert(nrows >= 0);
  assert(ncols >= 0);
  assert(ndata >= 0);
  sparse_t *A;
  if(!(A = malloc(sizeof(sparse_t))))
  {
    printf("Allocation failure in allocate_sparse - aborting.");
    exit(1);
  }

  A->m     = nrows;
  A->n     = ncols;
  A->nres  = ndata;
  if(!(A->ia = malloc((nrows+1)*sizeof(int))))
  {
    printf("Allocation failure in allocate_sparse - aborting.");
    exit(1);
  }
  /* Ensure that the size of a zero-dimensional matrix can be looked up. */
  A->ia[0] = 0;

  if(!(A->ja = malloc(ndata*sizeof(int))))
  {
    printf("Allocation failure in allocate_sparse - aborting.");
    exit(1);
  }

  if(!(A->a   = malloc(ndata*sizeof(double))))
  {
    printf("Allocation failure in allocate_sparse - aborting.");
    exit(1);
  }

  return A;
}


/*
 * Input arguments:
 *
 *     nrows    -
 *
 *     ncols    -
 *
 *     ndata    -
 *
 *     g        - Old sparse matrix
 *
 * Output:
 *
 *     return   - Sparse nrows x ncols matrix witgh room for ndata nonzeros.
 *
 */
/*---------------------------------------------------------------------------*/
sparse_t *sparse_realloc(sparse_t *g, int nrows, int ncols, int ndata)
/*---------------------------------------------------------------------------*/
{
  if (!g)
  {
    return sparse_alloc(nrows, ncols, ndata);
  }

  assert(nrows >= 0);
  assert(ncols >= 0);
  assert(ndata >= 0);

  g->m     = nrows;
  g->n     = ncols;
  g->nres  = ndata;
  if(!(g->ia = realloc(g->ia, (nrows+1)*sizeof(int))))
  {
    printf("Allocation failure in allocate_sparse - aborting.");
    exit(1);
  }
  g->ia[0] = 0;

  if(!(g->ja = realloc(g->ja, ndata*sizeof(int))))
  {
    printf("Allocation failure in allocate_sparse - aborting.");
    exit(1);
  }

  if(!(g->a   = realloc(g->a, ndata*sizeof(double))))
  {
    printf("Allocation failure in allocate_sparse - aborting.");
    exit(1);
  }

  return g;
}

/*
 * Input arguments:
 *
 *     to       - matrix pointer to write to, or pointer initialised with NULL
 *
 *     from     - sparse matric to read from.
 *
 * Output:
 *
 */
/*---------------------------------------------------------------------------*/
void sparse_copy (sparse_t **to, sparse_t *from)
/*---------------------------------------------------------------------------*/
{
  const int m = from->m;
  const int n = from->n;
  const int N = from->ia[m];
  if (!*to)
  {
    *to = sparse_alloc(m, n, N);
  }
  else if(((*to)->m != from->m) || ((*to)->ia[(*to)->m] < from->ia[from->m]))
  {
    *to = sparse_realloc(*to, m, n, N);
  }

  (*to)->m = from->m;
  (*to)->n = from->n;
  (*to)->nres = from->nres;

  memcpy( (*to)->ia, from->ia, (m+1)* sizeof(int)    );
  memcpy( (*to)->ja, from->ja,   N  * sizeof(int)    );
  memcpy( (*to)->a,  from->a,    N  * sizeof(double) );

}


/*
 * Input arguments:
 *
 *     g        - sparse matrix allocated with alloc_sparse
 *
 *     ncols    -
 *
 *     ndata    -
 *
 * Output:
 *
 *     return   - Sparse nrows x ncols matrix witgh room for ndata nonzeros.
 *
 */
/*---------------------------------------------------------------------------*/
void sparse_free (sparse_t *g)
/*---------------------------------------------------------------------------*/
{
  if (!g)
    return;

  if (g->ia)
  {
    free (g->ia);
    g->ia = NULL;
  }
  if (g->ja)
  {
    free (g->ja);
    g->ja = NULL;
  }
  if (g->a)
  {
    free (g->a);
    g->a = NULL;
  }
  free(g);
  return;
}



/*
 * Input arguments:
 *
 *     A        - pointer to sparse matrix
 *
 * Output:
 *
 *     A        - pointer to sparse matrix with transposed data.
 *                Internal storage may have changed due to reallocation of
 *                memory.
 *
 *     return   - A
 *
 */
/*---------------------------------------------------------------------------*/
sparse_t *sparse_transpose(sparse_t *A)
/*---------------------------------------------------------------------------*/
{

  int        i;
  int    k;
  int    N    = A->ia[A->m];
  int        rows = A->m;
  int        cols = A->n;

  int   *ia    = malloc ((cols+1) * sizeof (*ia ));
  int   *ja    = malloc (N        * sizeof (*ja ));
  double    *a     = malloc (N        * sizeof (*a  ));
  int       *count =  calloc(cols,      sizeof (*count));

  int pos;
#if 0
  for (i=0; i<rows; ++i)
  {
    for (k=A->ia[i]; k<A->ia[i+1]; ++k)
    {
      ++count [A->ja[k]];
    }
  }
#else
  int *ip = A->ja;
  for (k=A->ia[0]; k<A->ia[rows]; ++k)
  {
    ++count [*ip++];
  }

#endif

  ia[0] = 0;
  for (i=0; i<cols; ++i)
  {
    ia[i+1] = ia[i]+count[i];
  }

  memset(count, 0, cols*sizeof (*count));
  for (i=0; i<rows; ++i)
  {
    for (k=A->ia[i]; k<A->ia[i+1]; ++k, ++pos)
    {
      int c   = i;
      int r   = A->ja[k];

      pos     = ia[r]+count[r];
      ja[pos] = c;
      a [pos] = A->a[k];
      ++count [r];
    }
  }


  free (A->ia); A->ia = ia;
  free (A->ja); A->ja = ja;
  free (A->a ); A->a  = a;

  int tmp=A->m; A->m=A->n; A->n=tmp;
  free (count);

  return A;
}

/*---------------------------------------------------------------------------*/
sparse_t * sparse_zeros          (sparse_t *mat)
/*---------------------------------------------------------------------------*/
{
  int i;
  for (i=0; i<mat->nres; ++i)
  {
    mat->a [i] = 0.0;
  }
  return mat;
}


/*
 * Input arguments:
 *
 *     A        - Pointer to sparse matrix.
 *
 *     n        - number of (i,j,v)-triplets
 *
 *     row      - i
 *
 *     col      - j
 *
 *     value    - v
 *
 *     nrows    -
 *
 *     ncols    -
 *
 * Output:
 *
 *     return   - Sparse nrows x ncols matrix with (i,j,value) data
 *
 */
/*---------------------------------------------------------------------------*/
sparse_t * sparse_fill(sparse_t *A, const int n, int *row, int *col, double *value,
		       int nrows, int ncols)
/*---------------------------------------------------------------------------*/
{
  /*
    Init sparse matrix from list of (i,j,value)-triplets.  The dimension of the
    sparse matrix is the smallest that fits the data.
  */


  /* Find max entry in row and col */
  int i,j,p;
  int rows = 0;
  int cols = 0;
  for (i=0; i<n; ++i)
  {
    rows = rows > row[i] ? rows : row[i];
    cols = cols > col[i] ? cols : col[i];
  }

  /* Since we assume numbering from zero... */
  rows++;
  cols++;

  rows = rows > nrows ? rows : nrows;
  cols = cols > ncols ? cols : ncols;

  A = sparse_realloc (A, rows, cols, n);

  /* Sort first on row index, then on column index */
  sort (n, row, col, value);




  /* Count length of each row. */
  int *count = ( int * ) calloc (rows, sizeof(int));
  for (i=0; i<n; ++i)
  {
    ++count[row[i]];
  }

  /* fill rows */
  p         = 0;
  A->ia[0]  = 0;

  for (i=0; i<rows; ++i)
  {
    int this_row = row[p];
    for (j=0; j<count[i]; ++j)
    {
      assert ( row[p] == this_row );

      A->ja[p] = col  [p];
      A->a [p] = value[p];
      ++p;
    }

    /* Store position of next row */
    A->ia[i+1] = p;
  }

  free(count);

  return A;
}




/*
 * Input arguments:
 *
 *     A        - pointer to sparse matrix
 *
 *     f        - file pointer
 *
 * Output:
 *
 */
/*---------------------------------------------------------------------------*/
void sparse_display(sparse_t *A, FILE* f)
/*---------------------------------------------------------------------------*/
{

  int i;
  int k;

  double row[A->m];
  int    nz [A->m];
  fprintf(f, "\n");

  for (i=0; i<A->m; ++i)
  {
    for (k=0; k<A->m; ++k)
    {
      row[k] = 0;
      nz [k] = 0;
    }

    for (k=A->ia[i]; k<A->ia[i+1]; ++k)
    {
      row[A->ja[k]] = A->a[k];
      nz [A->ja[k]] = 1;
    }

    for (k=0; k<A->m; ++k)
    {
      if (nz[k])
	fprintf(f, "% 2.2f ", row[k]);
      else
	fprintf(f, "  .   ");
    }
    fprintf(f, "\n");
  }
}




/*
 * Input arguments:
 *
 *     A        -
 *
 *     fp       -
 *
 * Output:
 *
 *
 */
/*---------------------------------------------------------------------------*/
void sparse_write(sparse_t *A, FILE *fp)
/*---------------------------------------------------------------------------*/
{
  int i;
  int k;
  for(i=0; i<A->m; ++i)
  {
    for (k=A->ia[i]; k<A->ia[i+1]; ++k)
    {
      fprintf(fp, "%6d %6d %22.16lf\n", i, (int)A->ja[k], A->a[k]);
    }
  }
}




/*
 * Input arguments:
 *
 *     A        - pointer to quadratic sparse matrix
 *
 *     P        - permutation
 *
 * Output:
 *
 *     A        - row- and column permuted sparse matrix
 *
 */
/*---------------------------------------------------------------------------*/
void sparse_permute(sparse_t *A, int *P)
/*---------------------------------------------------------------------------*/
{
  /* Permute rows and columns of sparse matrix,
     such that A_{ij} = B_{P[i], P[j]}*/
  int    i,q;
  int k;
  int    m = A->m;
  if (m != A->n)
    {
      print("Error, cannot permute non-square matrix\n");
      exit(1);
    }

  sparse_t *B   = sparse_alloc(m, m, A->ia[m]);
  int      *invP = calloc(m, sizeof(int));

  /* Invert permutation */
  for (i=0; i<m; ++i)
    invP[P[i]] = i;

  /* Permuted row positions */
  B->ia[0]=0;
  for (i=0; i<m; ++i)
  {
    B->ia[i+1] = B->ia[i] + A->ia[P[i]+1]-A->ia[P[i]];
  }

  /* Permute columns */
  for (i=0; i<m; ++i)
  {
    const int pi = P[i];
    for (k=A->ia[pi], q=B->ia[i]; k<A->ia[pi+1]; ++k, ++q)
    {
      int j    = A->ja[k];
      B->ja[q] = invP[j];
      B->a [q] = A->a[k];
    }
  }

  /* Let A take over data structures of B */
  free (A->ia); A->ia = B->ia;
  free (A->ja); A->ja = B->ja;
  free (A->a ); A->a  = B->a;
  free (B);
  free (invP);
}





/*
 * Input arguments:
 *
 *     A        -
 *
 *     P        -
 *
 *     n_blocks -
 *
 * Output:
 *
 *     return   - 1 if A is block lower triangular
 *
 */
/*---------------------------------------------------------------------------*/
int sparse_check_block_triangularity(const sparse_t *A,       /* block-triangular?  */
			             const int      *P,       /* block pointers     */
			             const int      n_blocks)
/*---------------------------------------------------------------------------*/
{


  int i,j;
  int k;

  int is_block_triangular = 1;
  for (i=0; i<n_blocks; ++i)
  {
    /* block i: */
    for (j=P[i]; j<P[i+1]; ++j)
    {
      /* Check each member of the block */
      for (k=A->ia[j]; k<A->ia[j+1]; ++k)
      {
	/* No neighbour should be greater that the greatest in the same block. */
	if (((int)A->ja[k]) > P[i+1]-1)
	{
	  print("Strong component %d is not sorted correctly.\n", i);
	  is_block_triangular = 0;
	  break;
	}
      }
    }
  }

  return is_block_triangular;
}



/*
 * Input arguments:
 *
 *     a        -
 *
 *     b        -
 *
 * Output:
 *
 *     return   -
 *
 */
/*---------------------------------------------------------------------------*/
static int cmp(const void *a, const void *b)
/*---------------------------------------------------------------------------*/
{
  struct entry {
    int i;
    int j;
    double val;
  };

  const struct entry *A = (const struct entry *)a;
  const struct entry *B = (const struct entry *)b;

  if (A->i>B->i) return 1;
  else if(A->i == B->i  && A->j>B->j)       return 1;
  else return -1;

}




/*
 * Input arguments:
 *
 *     n        -
 *
 *     row      -
 *
 *     col      -
 *
 *     value    -
 *
 * Output:
 *
 *
 */
/*---------------------------------------------------------------------------*/
static void sort(const int n, int *row, int *col, double *value)
/*---------------------------------------------------------------------------*/
{

  typedef struct
  {
    int i;
    int j;
    double value;
  } element_t;

  element_t *e, *elms=malloc(n*sizeof(element_t));
  int i;
  for (i=0, e=elms; i<n; ++i, ++e)
  {
    e->i     = row  [i];
    e->j     = col  [i];
    e->value = value[i];
  }

  /* Sort elms */
  qsort((void*)elms, n, sizeof(element_t), &cmp);

  for (i=0, e=elms; i<n; ++i, ++e)
  {
    row  [i] = e->i;
    col  [i] = e->j;
    value[i] = e->value;
  }

  free(elms);
  return;
}






/*
 * Input arguments:
 *
 *     fp       -
 *
 *     S        -
 *
 * Output:
 *
 *
 */
/*---------------------------------------------------------------------------*/
int
sparse_read(FILE *fp, sparse_t **S)
/*---------------------------------------------------------------------------*/
{
  if (!fp)
  {
    print("In sparse_read:\n"
	  "Couldn't read filepointer\n");
    return 1;
  }

  int    bufsz  = 1000;
  int    *row   = (int *)   malloc( bufsz * sizeof(int)   );
  int    *col   = (int *)   malloc( bufsz * sizeof(int)   );
  double *value = (double *)malloc( bufsz * sizeof(double));

  /* Read data from file */
  int m, n;
  fscanf(fp, "%d %d\n", &m, &n);

  int pos = 0;
  while (fscanf(fp, "%d %d %lf\n", row+pos,col+pos,value+pos) == 3)
  {
    /* If data is numbered from 1 */
#if 0
    --row[pos];
    --col[pos];
#endif
    ++pos;

    /* Increase buffers? */
    if (pos >= bufsz)
    {
      bufsz = (bufsz *1.5);/* lround(bufsz *1.5); */
      row   = (int *)    realloc( row,   bufsz * sizeof(int)   );
      col   = (int *)    realloc( col,   bufsz * sizeof(int)   );
      value = (double *) realloc( value, bufsz * sizeof(double));
    }
  }

  /* Resize storage to fit data. */
  bufsz = pos;
  row   = (int *)    realloc( row,   bufsz * sizeof(int)   );
  col   = (int *)    realloc( col,   bufsz * sizeof(int)   );
  value = (double *) realloc( value, bufsz * sizeof(double));

  /* Make sparse matrix from unsorted data */
  *S = NULL;
  *S = sparse_fill(*S, bufsz, row, col, value, m, n);
  (*S)->m = m;
  (*S)->n = n;
  assert ( (*S)->m == (*S)->n );

  /* Free storage */
  free (row);
  free (col);
  free (value);
  return 0;
}




/*
 * Input arguments:
 *
 *     A        -
 *
 * Output:
 *
 *     return   - 1 if A is strictly lower triangular
 *
 */
/*---------------------------------------------------------------------------*/
int
sparse_is_triangular(sparse_t *A)
/*---------------------------------------------------------------------------*/
{
  int n = A->m;
  int i,k;

  for (i=0; i<n; ++i)
  {
    for (k=A->ia[i]; k<A->ia[i+1]; ++k)
    {
      int j = A->ja[k];

      if (j>i)
	return 0;
    }
  }
  return 1;
}




/*
 * Input arguments:
 *
 *     A        -
 *
 *     b        -
 *
 *     x        -
 *
 * Output:
 *
 *
 */
/*---------------------------------------------------------------------------*/
void
sparse_solve_triangular(sparse_t *A, double *b, double *x)
/*---------------------------------------------------------------------------*/
{
  assert (A->m == A->n);
  int n = A->m;
  int i,k;
  for (i=0; i<n; ++i)
  {
    double ajxj  = 0.0;
    double aii   = 0.0;
    for (k=A->ia[i]; k<A->ia[i+1]; ++k)
    {
      int j = A->ja[k];

      if (j<i)
	ajxj += A->a[k]*x[j];
      else /* i==j */
	aii = A->a[k];
    }

    if (fabs(aii)<EPS)
      print("In sparse_solve_lower_triangular: Matrix is singular.\n");


    x[i] = (b[i] - ajxj)/aii;
  }
}


/*---------------------------------------------------------------------------*/
static void
solve_full_diag_block(int begin, int end, sparse_t *A, double *B, double *x)
/*---------------------------------------------------------------------------*/
{
  int n = end-begin;
  double *a = calloc(n*n,  sizeof(double));
  double *b = malloc(n*    sizeof(double));
  int i,j,k;

  for (i=begin; i<end; ++i)
    b [i-begin] = B [i];


  for (i=begin; i<end; ++i)
  {
    int II = i-begin;
    for (k=A->ia[i]; k<A->ia[i+1]; ++k)
    {

      j = A->ja[k]-begin;	  assert(j<end);
      if (j<II)
	b[II] -= A->a[k]*x[begin+j];
      else /* i==j */
	a[j*n+II] += A->a[k];
    }
  }

  lu(n, a, b);
  memcpy(x+begin, b, n*sizeof *b);

  free(b);
  free(a);
}



/*---------------------------------------------------------------------------*/
static void
solve_sparse_diag_block(int begin, int end, sparse_t *A, double *B, double *x)
/*---------------------------------------------------------------------------*/
{
  int       n = end-begin;
  int       i,j,k;
  int       nnz = A->ia[end]-A->ia[begin];
  sparse_t *a = sparse_alloc(n, n, nnz);
  double   *b = calloc(n, sizeof *b);
  int       pos = 0;

  for (i=begin; i<end; ++i)
    b [i-begin] = B [i];

  a->ia[0] = 0;
  for (i=begin; i<end; ++i)
  {
    int I = i-begin;
    for (k=A->ia[i]; k<A->ia[i+1]; ++k)
    {

      j = A->ja[k]-begin;	  assert(j<end);
      if (j<0)
	b[I] -= A->a[k]*x[begin+j];
      else /* i==j */
      {
	a->ja[pos]  = j;
	a->a [pos]  = A->a[k];
	++pos;
      }
    }
    a->ia[i+1-begin]=pos;
  }
  /* sparse_display(a, stderr);fprintf(stderr, "\n"); */

  sparse_solve(a, b, x+begin);

  free(b);
  sparse_free(a);
}

/*---------------------------------------------------------------------------*/
static void
solve_scalar_diag_block(int begin, int end, sparse_t *A, double *B, double *x)
/*---------------------------------------------------------------------------*/
{
  int size = end-begin;
  int k;
  assert(size == 1);

  double ajxj  = 0.0;
  double aii   = 0.0;
  for (k=A->ia[begin]; k<A->ia[end]; ++k)
  {
    int j = A->ja[k];

    if (j<begin)
    {
      ajxj += A->a[k]*x[j];
    }
    else /* i==j */
    {
      aii = A->a[k];
    }
  }

  if (fabs(aii)<EPS)
  {
    print("In sparse_solve_lower_triangular: Matrix is singular.\n");
  }

  x[begin] = (B[begin] - ajxj)/aii;
}

/*
 * Input arguments:
 *
 *     A        -
 *
 *     b        -
 *
 *     x        -
 *
 *     nc       - number of strong components
 *
 *     C        - Int vector of pointers to strong components
 *
 * Output:
 *
 *
 */
/*---------------------------------------------------------------------------*/
void sparse_solve_block_triangular(sparse_t *A, double *b, double *x,
				   int nc, int *C)
/*---------------------------------------------------------------------------*/
{
  assert (A->m == A->n);
  int i;
#if 0
  int begin;
  int end;

  i   = 0;
  end = C[i++];
  while (i<nc)
  {
    begin = end;
    end = C[i];

    if (end-begin == 1)
    {
      /* Pass as many 1x1 diagonal blocks as possible to scalar solver*/
      while (C[i]-C[i-1] == 1)
      {
	++i;
      }
      end = C[i];
      solve_scalar_diag_block(begin, end, A, b, x);
    }
    else if (end-begin < 100)
    {
      solve_full_diag_block   (begin, end, A, b, x);
    }
    {
      solve_sparse_diag_block (begin, end, A, b, x);
    }


  }
#else
  /* If A is simply lower triangular, the following code has less overhead We
   * could maybe combine the two by locating the loops first, and use a less
   * overhead on the long sequences of simple 1x1 blocks
   */
  if (sparse_is_triangular(A))
  {
    sparse_solve_triangular(A, b, x);
    return;
  }

  for (i=0; i<nc; ++i)
  {
    int size = C[i+1]-C[i];

    if (size == 1)
    {
      solve_scalar_diag_block(C[i],C[i+1], A, b, x);
    }
    else if (size < 100)
    {
      solve_full_diag_block   (C[i], C[i+1], A, b, x);
    }
    {
      solve_sparse_diag_block (C[i], C[i+1], A, b, x);
    }


#endif


  }
}





/*
 * Input arguments:
 *
 *     A        -
 *
 *     b        -
 *
 *     x        -
 *
 * Output:
 *
 *
 */
/*---------------------------------------------------------------------------*/
void
sparse_solve(sparse_t *A, double *b, double *x)
/*---------------------------------------------------------------------------*/
{
  assert (A->m == A->n);

  if (sparse_is_triangular(A))
  {
    sparse_solve_triangular(A, b, x);
    return;
  }


  int n = A->m;

  /* UMFPACK expect compressed-sparse-column format. */
  sparse_t *B = NULL;  sparse_copy(&B, A);
  B = sparse_transpose(B);

  double *null = (double *) NULL ;

  void *Symbolic, *Numeric ;

  umfpack_di_symbolic (n, n, B->ia, B->ja, B->a, &Symbolic, null, null) ;
  umfpack_di_numeric  (B->ia, B->ja, B->a, Symbolic, &Numeric, null, null) ;
  umfpack_di_free_symbolic (&Symbolic);


  umfpack_di_solve (UMFPACK_A, B->ia, B->ja, B->a, x, b, Numeric, null, null) ;
  umfpack_di_free_numeric (&Numeric);

  sparse_free(B);
}




/*---------------------------------------------------------------------------*/
void
sparse_to_full(sparse_t *A, double *F, int by_rows)
/*---------------------------------------------------------------------------*/
{
  const int m = A->m;
  const int n = A->n;

  int i,j,k;

  if (by_rows) /* column index runs fastest */
  {
    for (i=0; i<A->m; ++i)
    {
      for (j=0; j<n; ++j)
      {
	F[j]=0.0;
      }

      for (k=A->ia[i]; k<A->ia[i+1]; ++k)
      {

	j = A->ja[k];
	F[j] += A->a[k];
      }
      F = F+n;
    }
  }
  else /* row index runs fastest */
  {
    for (k=0; k<m*n; ++k)
    {
      F[k]=0.0;
    }
    for (i=0; i<A->m; ++i)
    {
      for (k=A->ia[i]; k<A->ia[i+1]; ++k)
      {

	j = A->ja[k];
	F[j*m+i] += A->a[k];
      }
    }
  }
}
