/*------------------------------------------------------------
 File       : sparse.h
 Created    : Thu Aug  7 16:10:06 2008
 Author     : Jostein R. Natvig <Jostein.R.Natvig@sintef.no>
 Revision   : $Id: sparse.h 136 2008-08-26 13:07:00Z jrn $


 Description


 Changes

------------------------------------------------------------*/

#ifndef _SPARSE_H_
#define _SPARSE_H_

typedef struct
{
  int     m;        /* Number of rows                                 */
  int     n;        /* Number of columns                              */
  int     nres;     /* Reserved space for non-zeros                   */
  int    *ia;       /* Position of row i                              */
  int    *ja;       /* Neighbour vertices (column number, j)          */
  double * a;       /* Flux to neighbour vertices (value in (i,j))    */
} sparse_t;

sparse_t * sparse_alloc    (int nrows, int ncols, int ndata);
sparse_t * sparse_realloc  (sparse_t *g, int nrows, int ncols, int ndata);
void       sparse_free     (sparse_t *g);
void       sparse_copy     (sparse_t **to, sparse_t *from);

sparse_t * sparse_transpose(sparse_t *A);
sparse_t * sparse_zeros    (sparse_t *A);

sparse_t * sparse_fill     (sparse_t *A, const int n, int *row, int *col,
			    double *value, int nrows, int ncols);
void       sparse_display  (sparse_t *A, FILE* f);
void       sparse_permute  (sparse_t *A, int *rowP);

int        sparse_is_triangular(sparse_t *A);
int        sparse_check_block_triangularity(const sparse_t *A, const int *P,
					    const int n_blocks);
void       sparse_write    (sparse_t *A, FILE *fp);
int        sparse_read     (FILE *fp, sparse_t **S);
void       sparse_to_full  (sparse_t *A, double *F, int by_rows);

void       sparse_solve    (sparse_t *A, double *b, double *x);
void       sparse_solve_triangular       (sparse_t *A, double *b, double *x);
void       sparse_solve_block_triangular (sparse_t *A, double *b, double *x,
					  int nc, int *C);

#endif
