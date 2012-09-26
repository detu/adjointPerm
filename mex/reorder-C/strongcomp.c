/*------------------------------------------------------------
 File       : strongcomp.c
 Created    : Wed Aug  6 10:42:05 2008
 Author     : Jostein R. Natvig <Jostein.R.Natvig@sintef.no>
 Revision   : $Id$


 Description
 Tarjans algorithm for finding strong components and topological
 sequence of a directed graph.  The graph is encoded in compressed
 sparse row format.

 There is a weakness in this implementation.  For a vertex with a
 very large number n of descendants results in an n^2 complexity.

 Changes

------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "sparse.h"
#include "strongcomp.h"


#if MEXFILE

#include <mex.h>
#define print mexPrintf
extern int interrupted;

#else

#define print printf
#define interrupted

#endif

#define TEST 1



static int
min(int a, int b){ return a<b? a: b;}


/*
 * Input arguments:
 *
 *     A        - n-vector to be permuted.
 *
 *     work     - n-vector for temporary storage.
 *
 *     p        - n-vector holding permutation such that
 *                vec[i] becomes vec[p[i]].
 *
 * Output:
 *
 *     rowP     - n-vector of enpoints of each strong component.
 *
 *     P        - n-vector holding topological ordering.
 *
 *     return   - number of strong components.
 */

/*---------------------------------------------------------------------------*/
int reorder(sparse_t *A, int* rowP, int *P, void *work)
/*---------------------------------------------------------------------------*/
{
  assert(A->m == A->n);
  int ncomp;

  strongcomp(A->m, A->ia, A->ja, rowP, P, &ncomp, work);

  return ncomp;
}


/*
 * Input arguments:
 *
 *     nc       - number of components
 *
 *     comp     - pointer to start of components, i.e., component
 *                number i has nodes comp[i], ..., comp[i+1]-1.
 *
 *     P        - Node permutation. If P!=NULL, Components become
 *                P[comp[i]], ..., P[comp[i]-1].
 *
 *     func     - function f(begin, end, data)
 *
 *     data     - void pointer to some data
 *
 * Output:
 *
 *     return   - 0 if success, nonzero  if error.
 *
 * NOTE:
 *
 *     We currently assume that the system is actually permuted to triangular shape
 */
#include "system.h"
/*---------------------------------------------------------------------------*/
int
foreach_component(int   nc,                    /* num components            */
		  int  *comp,                  /* begin and end of comp     */
		  int  *order,                 /* sequence                  */
		  int (*solvescalar)(int, int, void *, void*), /* one cell  */
		  int (*solvesystem)(int, int, void *, void*), /* may cells */
		  void *data,
		  void *report)
/*---------------------------------------------------------------------------*/
{
  /* For each component, apply func to data */
  int i;

  if (order==NULL)
  {
    for (i=0; i<nc; ++i)
    {
      int begin  = comp[i];
      int end    = comp[i+1];
      int err;

      if ((end-begin) == 1) err = solvescalar(begin, end, data, report);
      else                  err = solvesystem(begin, end, data, report);

      if (err)
      {
	print("Error occured in component %d (%d...%d).\n",
	      i, begin, end-1);

	return 1;
      }


#if MEXFILE
      if (interrupted)
      {
	print("interrupted");
	return 1;
      }
#endif
    }
  }
  else
  {
    print("order!=NULL not currently supported in 'foreach_component'");
    return 1;
  }
  return 0;
}




/*---------------------------------------------------------------------------*/
int write_example_to_file(sparse_t *A, int* rowP)
/*---------------------------------------------------------------------------*/
{
  extern int interrupted;
  /* Write matrix and sequence to file*/
  char matrix_file  [16] = "A.txt";
  char sequence_file[16] = "A.sequence";
  FILE *fp;
  if (!(fp = fopen(matrix_file, "w")))
    {
      print("Error! Could not open file %s", matrix_file);
      return 1;
    }

  int i,k;
  for (i=0; i<A->m; ++i)
  {
    for (k=A->ia[i]; k<A->ia[i+1]; ++k)
    {
      const int    j = A->ja[k];
      const double a = A->a [k];
      fprintf(fp, "%d %d %f\n", i, j, a);
    }
  }
  fclose(fp);

  if (!(fp = fopen(sequence_file, "w")))
    {
      fprintf(stderr, "Error! Could not open file %s", matrix_file);
      exit(1);
    }
  for (k=0; k<A->m; ++k)
    fprintf(fp, "%d\n", rowP[k]);

  fclose(fp);

  return 0;
}


/*
 * Input arguments:
 *
 *     size     - Number of nodes.
 *
 *     ia       - Row pointers for sparse adjacency matrix.
 *
 *     ja       - Column numbers for sparse adjacency matrix.
 *
 *     P        - Int vector of length SIZE to hold topological sequence.  Used
 *                temporary as stack.
 *     Q        - Int vector of length SIZE+1 to hold start end end pointers of
 *                each strong component.
 *     ncomp    -
 *
 *     work     - Poiniter to block of memory size*(2*sizeof(int)+sizeof(char))
 *                long used as temporary storage.
 * Output:
 *
 *     P        - Topological sequence.
 *
 *     Q        - Pointers to beginning (and end) of each strong component.
 *
 *     ncomp    - Number of strong components
 */



/* Improved (or obfuscated) version uses less memory
 *
 * Use end of P and Q as stack:
 *   push operation is *s--=elm,
 *   pop  operation is elm=*++s,
 *   peek operation is *(s+1)
 */

void
strongcomp(int size, int *ia, int *ja, int *P, int *Q, int *ncomp, int *work)
/*---------------------------------------------------------------------------*/
{
  /* To eliminate malloc/free */
  int   *time  = work;
  int   *link  = (int *)  time + size;
  char  *inssc = (char*) (link + size); /* byte... */

  /* Stacks */
  int *s      = Q + size;      /* visit stack            */
  int *ssc    = P + size-1;    /* strong component stack */

  memset(time, 0, size * (sizeof(*time) + sizeof(*link) + sizeof(*inssc)));
  memset(P,    0, size * sizeof *P   );
  memset(Q,    0, size * sizeof *Q   );



  int  thetime = 0;
  int  count   = 0;
  int *bottom  = s;
  int  seed;

  *ncomp = 0;
  *Q++   = 0;
  for (seed=0; seed<size; ++seed)
  {
    /* Skip if seed already visited */
    if (time[seed]>0) continue;

    /* Push seed on stacks s and ssc */
    *s--        = seed;
    *ssc--      = seed;

    inssc[seed] = 1;

    time [seed] = ++thetime;
    link [seed] =   thetime;

    while (s != bottom)
    {
      /* Get vertex at top of stack s */
      int c  = *(s+1);
      int p;
      for (p=ia[c]; p<ia[c+1]; ++p)
      {
	int n = ja[p];
	if(time[n] == 0)
	{
	  /* Push vertex n on stacks s and ssc */

	  assert (ssc >= P);
	  assert (s   >= Q);
	  *s--     = n;
	  *ssc--   = n;

	  assert(inssc[n] == 0);
	  inssc[n] = 1;

	  time [n] = ++thetime;
	  link [n] =   thetime;

	  link [c] = min(link[c], link[n]);
	  break;
	}

	if (inssc[n])
	{
	  link[c] = min(link[c], link[n]);
	}
      }

      /* All descendants has been visited */
      if ( p == ia[c+1] )
      {
	/* Pop vertex from visit stack */
	int tmp = *(++s);
	assert(tmp == c);


	/* If c is root of subtree that is a strong component */
	if (link[c] == time[c])
	{

	  int v;
	  /* Pop vertices belongin to same strong-component from ssc stack */
	  do{
	    v        = *++ssc;
	    inssc[v] = 0;

	    assert (( ssc >= P  )  |   /* There's still room on the ssc stack */
		    ( s == bottom ) );/* or we're finished */

	    *P++     = v;
	    ++count;
	  }while (v != c);

	  /* Store end of strong component */
	  assert (s >= Q);
	  *Q++ = count;
	  ++*ncomp;

	  /* Loop removal could be implemented here: untangle each completely
	     explored strong component */

	  /* deloop(count, ia, ja, P, Q, work); */

	}
	else
	{
	  /* c belongs to a strong component but is not the root. */
	}
      }
    }
  }
}




#if TEST
#include "sparse.h"


/* Test code.  Reads a sparse matrix from a file specified on the   */
/* command line with file format "m n\n i1 j1 v1\n i2 j2 v2\n ...". */
/* Computes reorder sequence                                        */
int main (int argc, char *argv [])
{
  sparse_t *S;

  FILE *fp;

  if ((fp = fopen(argv[1],"r")))
  {
    sparse_read(fp, &S);
    fclose(fp);
  }

  /* sparse_display(S, stderr); */
  int *P = calloc(S->m, sizeof(int));
  int *Q = calloc(S->m, sizeof(int));
  int *work  = malloc (S->m* (2*sizeof *work + sizeof (char)));

  int ncomp;
  strongcomp(S->m, S->ia, S->ja, P, Q, &ncomp, work);



  int i,j;
  for (i=0; i<ncomp; ++i)
  {
    for (j=Q[i]; j<Q[i+1]; ++j)
    {
      fprintf(stderr, "%d\n", P[j]);
    }
  }


  free(P);
  free(Q);
  sparse_free(S);

  return 0;
}
#endif
