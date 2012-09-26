#include <stdlib.h>
#include <stdio.h>
#include "lu.h"

void dgetrf_ (int *M,int *N,double *A,int *LDA,
	      int *IPIV, int *INFO);

void dgetrs_ (char *trans, int *N, int *NRHS, double *A,int *LDA,
	      int *IPIV, double *B, int *LDB, int *INFO);
void dgesv_ (int*, int*, double*, int*, int*, double*, int*, int*);
/******************************************************************************
 *                                                                            *
 * LAPACK wrappers                                                            *
 *                                                                            *
 *****************************************************************************/
/* LAPACK linear solver */
/*---------------------------------------------------------------------------*/
int
lu (int N, double *A, double *B)
/*---------------------------------------------------------------------------*/
{
  int info;
  int lda  = N;
  int ldb  = N;
  int NRHS = 1;
  int ipiv [N];

  /* Simple LAPACK LU-routine */
  dgesv_ (&N, &NRHS, A, &lda, ipiv, B, &ldb, &info);

  if (info > 0){
    fprintf(stdout, "\nNonzero return value from LAPACK: %d\n", info);
    /* print_matrix(stdout, N, N, A); */
  }
  return (info > 0);
}


/* LAPACK LU factorisation */
/*---------------------------------------------------------------------------*/
int
lufac (int N, double *A)
/*---------------------------------------------------------------------------*/
{
  int info;
  int lda  = N;
  int ipiv [N];
  dgetrf_ (&N, &N, A, &lda, ipiv, &info);
  return info;
}


/* LAPACK Solve using LU factorisation */
/*---------------------------------------------------------------------------*/
int
lusolve (int N, double *A, double *B)
/*---------------------------------------------------------------------------*/
{
  char trans = 'N';
  int info;
  int lda  = N;
  int ldb  = N;
  int NRHS = 1;
  int ipiv [N];
  dgetrs_ (&trans, &N, &NRHS, A, &lda, ipiv, B, &ldb, &info);
  return info;
}


/*---------------------------------------------------------------------------*/
void
print_matrix(int m, double **A)
/*---------------------------------------------------------------------------*/
{
  int i,j;
  printf("\n");
  for (i = 0; i < m; i++)
  {
    printf("[");
    for (j = 0; j < m; j++)
      printf("% 4.4g ", A[i][j]);
    printf("]\n");
  }
}
/*---------------------------------------------------------------------------*/
void
transpose_matrix(int m, double **A)
/*---------------------------------------------------------------------------*/
{
  int i,j;
  for (i = 0; i < m; i++)
  {
    for (j = 0; j < i; j++)
    {
      double tmp = A[i][j];
      A[i][j] = A[j][i];
      A[j][i] = tmp;
    }
  }
}


/*---------------------------------------------------------------------------*/
void
print_vector(int m, double *b)
/*---------------------------------------------------------------------------*/
{
  int j;
  printf("\n");
  printf("[");
  for (j = 0; j < m; j++)
    printf("%f ", b[j]);
  printf("]\n");

}


#if 0
/*---------------------------------------------------------------------------*/
int main()
/*---------------------------------------------------------------------------*/
{
  double **A,*b;
  int *ipiv;
  const int m = 3;

  NEW_MATRIX(A, double, m, m);
  NEW_ARRAY (b, double, m);
  NEW_ARRAY (ipiv, int, m);

  int i,j;
  for (i = 0; i < m; i++)
    for (j = 0; j < m; j++)
    {
      A[i][j]=i+j;
      b[i] = i;
    }

  A[2][2] = 5;
  print_matrix(m,A);
  /* Simple driver. */
  dgetrf(m, m, A, ipiv);
  print_matrix(m,A);

  print_vector(m, b);

  dgetrs(m, A, ipiv, b);

  print_vector(m, b);
  print_matrix(m,A);



}
#endif
