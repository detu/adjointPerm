/*------------------------------------------------------------
 File       : ispu.c
 Created    : Thu Aug  7 13:37:53 2008
 Author     : Jostein R. Natvig <Jostein.R.Natvig@sintef.no>
 Revision   : $Id: ispu.c 192 2008-10-30 07:13:50Z jrn $


 Description
 Mex interface to reordered first-order upwind solver with(out)
 compressibility.  Version 2.0

 Changes

------------------------------------------------------------*/
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#define ENABLE_MEXERRMSGTXT 1
#define CHECK_DIVERGENCE    0

#include <string.h>
#include <mex.h>
#include <matrix.h>
#include <math.h>
#include <stdio.h>
#include <stdarg.h>

#include "sparse.h"
#include "strongcomp.h"
#include "system.h"
#include "utils.h"
#include "mexutils.h"

#include "today.h"

extern int    interrupted;
extern double mobratio;



/* ---- Settings ----------------------------- */
struct options
{
  double mobility_ratio;
  int    max_iterations;
  int    substeps;
  double tolerance;
  int    verbosity;
  char   *scalar_solver;
};

enum field {MOBILITY_RATIO,
	    MAX_ITERATIONS,
	    TOLERANCE,
	    VERBOSITY,
	    SCALAR_SOLVER};



/* ----- Statistics and diagnostic output ------*/
static const int strlength = 16384;
struct report
{
  double average_iterations;
  int    max_iterations;
  int    num_zeroiterations;
  int    max_loopsize;
  int    num_loops;
  int    num_looped_cells;
  double internal_time;
  char   message[16384];
};



static int  isFull          (const mxArray *);
static int  insaneInput     (const int, const mxArray**);
static void get_options     (struct options *, int[], const mxArray *);
static int  solve           (double, sparse_t*, double*, double*, double*,
			     struct options*, struct report*);
#if CHECK_DIVERGENCE
static int  checkDivergence (const sparse_t*, const double*);
#endif





/*---------------------------------------------------------------------------*/
void  export_report(const struct report *in, mxArray **out)
/*---------------------------------------------------------------------------*/
{

  const mwSize  nfields       = 8;
  const char   *fieldnames[8] = {"average_iterations",
				 "max_iterations",
				 "num_zeroiterations",
				 "max_loopsize",
				 "num_loops",
				 "num_looped_cells",
				 "internal_time",
				 "message"};
  mxArray *arr;
  *out = mxCreateStructMatrix(1, 1, nfields, fieldnames);


  arr = exportArray(&in->average_iterations, sizeof(double),
		    1, 1, mxDOUBLE_CLASS);
  mxSetField (*out, 0, "average_iterations", arr);



  arr = exportArray(&in->max_iterations,     sizeof(int),
		    1, 1, mxINT32_CLASS);
  mxSetField (*out, 0, "max_iterations",     arr);


  arr = exportArray(&in->num_zeroiterations, sizeof(int),
		    1, 1, mxINT32_CLASS);
  mxSetField (*out, 0, "num_zeroiterations",     arr);



  arr = exportArray(&in->max_loopsize,       sizeof(int),
		    1, 1, mxINT32_CLASS);
  mxSetField (*out, 0, "max_loopsize",       arr);



  arr = exportArray(&in->num_loops,          sizeof(int),
		    1, 1, mxINT32_CLASS);
  mxSetField (*out, 0, "num_loops",          arr);



  arr = exportArray(&in->num_looped_cells,   sizeof(int),
		    1, 1, mxINT32_CLASS);
  mxSetField (*out, 0, "num_looped_cells",   arr);



  arr = exportArray(&in->internal_time,      sizeof(double),
		    1, 1, mxDOUBLE_CLASS);
  mxSetField (*out, 0, "internal_time",      arr);



  arr = exportArray(in->message,             sizeof(char),
		    1, strlen(in->message), mxCHAR_CLASS);
  mxSetField (*out, 0, "message",            arr);



}



/*---------------------------------------------------------------------------*/
void mexFunction(int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
/*---------------------------------------------------------------------------*/
{
  if (insaneInput(nrhs, prhs)) return;





  /* Install handler for SIGINT signal */
  sighandler_t oldInterruptHandler = signal(SIGINT, interruptHandler);
  interrupted    = 0;



  /* Fetch input and allocate output */
  double   tf;      copyMatlabScalar(&tf,      prhs[3]);


  sparse_t *flux  = getSparse        (prhs[0]);
  double   *pv    = copyMatlabVector (prhs[1]);
  plhs[0]                = mxDuplicateArray (prhs[2]);
  double   *sat   = mxGetPr          (plhs[0]);
  double   *q     = copyMatlabVector (prhs[4]);


  struct options opt = {0};
  int fieldsize[5];


  get_options(&opt, fieldsize, prhs[5]);

/*   mexPrintf("Solve using %s\n", opt.scalar_solver); */


  struct report r;
  solve(tf, flux, pv, q, sat, &opt, &r);

  /* Free local copies of Matlab data */
  free(q);
  free(pv);
  sparse_free(flux);




  if (nlhs == 2)
  {
    export_report (&r, &plhs[1]);
  }


  /* Restore old SIGINT handler */
  signal(SIGINT, oldInterruptHandler);

#if ENABLE_MEXERRMSGTXT
  /* Force an interrupt in Matlab */
  if ( interrupted )
    {
      mexErrMsgTxt("User interrupt");
    }
#endif
}

/*---------------------------------------------------------------------------*/
void get_options (struct options *out, int num_elements[5], const mxArray *in)
/*---------------------------------------------------------------------------*/
{
  const int   n = 6;
  int         offset[6];
  char       *names [6] = {"mobility_ratio",
			   "max_iterations",
			   "substeps",
			   "tolerance",
			   "verbosity",
			   "scalar_solver"};
  enum returnType rtype[6] = {VALUE,
			      VALUE,
			      VALUE,
			      VALUE,
			      VALUE,
			      POINTER};


  offset[0] = offsetof(struct options, mobility_ratio);
  offset[1] = offsetof(struct options, max_iterations);
  offset[2] = offsetof(struct options, substeps);
  offset[3] = offsetof(struct options, tolerance);
  offset[4] = offsetof(struct options, verbosity);
  offset[5] = offsetof(struct options, scalar_solver);


  /*   Use this instead */
  mxClassID  id[6]= {mxDOUBLE_CLASS,
		     mxINT32_CLASS,
		     mxINT32_CLASS,
		     mxDOUBLE_CLASS,
		     mxINT32_CLASS,
		     mxCHAR_CLASS};


  import_struct(out, n, num_elements, names, offset, id, rtype, in);
}

/* void free_options(struct options opt); */

/*---------------------------------------------------------------------------*/
void loop_statistics(int *P, int nc, int *num_looped_cells,
		     int *num_loops, int *max_loopsize)
/*---------------------------------------------------------------------------*/
{
  int i;
  *num_looped_cells = 0;
  *num_loops        = 0;
  *max_loopsize     = 0;


  for (i=0; i<nc; ++i)
  {
    if (P[i+1]-P[i] > 1)
    {
      ++(*num_loops);
      int loopsz   = P[i+1]-P[i];
      *num_looped_cells  += loopsz;
      *max_loopsize = *max_loopsize > loopsz ? *max_loopsize : loopsz;
    }
  }

}


/*---------------------------------------------------------------------------*/
int Message (char **buf, int *bufsize, const char *format, ...)
/*---------------------------------------------------------------------------*/
{
  int size;

  va_list ap;
  va_start(ap, format);

  size = vsprintf(*buf, format, ap);
  mexPrintf("%s", *buf);
  *buf     += size;
  *bufsize += size;
  va_end(ap);

  return size;
}
/*
 * Input arguments:
 *
 *     tf       - final time
 *
 *     dt       - magnitude of time step
 *
 *     flux     - nxn sparse matrix of (positive) in-fluxes.
 *                (Confusing: The same as Fin in computeFlux)
 *
 *     volume   - n-vector holding pore volumes
 *
 *     sources  - n-vector holding sources and sinks for each cell
 *
 *     u        - n-vector holding current satruations S.
 *
 *     PT       - n-vector holding dP/dt*ci*phi
 *
 *     PG       - n-vector holding c_i*v_i*\grad P
 *
 * Output:
 *
 *     u        - n-vector of new saturations S.
 *
 *
 * The implicit upwind scheme for incompressible twopahse transport
 *
 *    s^n - s^(n-1) + PT*s^n + dt/volume*flux*f(s)+PG*f(s) = sources
 *
 * is used to advance the solution from t=0 to t=tf in time steps of
 * size dt.
 */
/*---------------------------------------------------------------------------*/
static int
solve(double          tf,
      sparse_t       *flux,
      double         *volume,
      double         *source,
      double         *u,
      struct options *opt,
      struct report  *r)
/*---------------------------------------------------------------------------*/
{
  double t;


  mobratio  = opt->mobility_ratio;
  char *str = r->message;
  int  len  = strlength;

  /* --- banner --- */
  Message(&str, &len, "\n");
  Message(&str, &len, "------------------------------------------------------------\n");
  Message(&str, &len, "Implicit upwind solver for incompressible two-phase flow  \n");
  Message(&str, &len, "%s,    maintainer jrn  \n", today);
  Message(&str, &len, "------------------------------------------------------------\n");




  const int n_cells = flux->m;
  int      *rowP    = malloc (  n_cells    * sizeof(int));
  int      *P       = malloc ( (n_cells+1) * sizeof(int));
  double   *work    = malloc (  5*n_cells  * sizeof(double));

  /* Permute flux */
  Message (&str, &len, "Reorder               :\t");
  tic(0);
  int nc = reorder(flux, rowP, P, work);
  Message (&str, &len, "Time elapsed: %10.7f seconds\n", toc(0));

  Message (&str, &len, "Permute flux          :\t");
  tic(0);
  sparse_permute(flux, rowP);
  Message (&str, &len, "Time elapsed: %10.7f seconds\n", toc(0));
  Message (&str, &len, "\n");

  Message (&str, &len, "Permute pv, q, s      :\t");
  tic(0);

  permute_dvector(volume, work, rowP, n_cells);
  permute_dvector(source, work, rowP, n_cells);
  permute_dvector(u,      work, rowP, n_cells);
  Message(&str, &len, "Time elapsed: %10.7f seconds\n", toc(0));



  if (!sparse_check_block_triangularity(flux, P, nc))
  {
    mexWarnMsgTxt ("The permuted matrix of fluxes is not block triangular\n");
  }


  /* Write report */
  Message (&str, &len, "\n"); fflush(stdout);
  Message (&str, &len, "Number of fluxes      :    %8d\n",   (int)flux->ia[flux->m]);
  Message (&str, &len, "Number of cells       :    %8d\n",   flux->m);
  Message (&str, &len, "Subdomains            :    %8d\n\n", nc);



  t  = 0.0;
  /* dt = dt < tf - t ? dt : tf - t; */



  system_t sys = {0};
  init_system(&sys, tf/opt->substeps, flux, source, volume, u,
	      opt->mobility_ratio, opt->max_iterations, opt->tolerance, work);

  /* Init report */
  report_t rep = {0};
  rep.iterations = calloc(flux->m, sizeof(int));
  rep.type       = calloc(flux->m, sizeof(char));


  Message (&str, &len, "Mobility ratio        :       %2.2f\n", opt->mobility_ratio );
  Message (&str, &len, "Max iterations        :    %8d\n",      opt->max_iterations );
  Message (&str, &len, "Tolerance             :       %2.2e\n", opt->tolerance);
  Message (&str, &len, "Computing solution");
  Message (&str, &len, "\n");




  r->internal_time = 0.0;
  while ((t < tf) && !interrupted)
  {
    Message (&str, &len, "Time step (%2.2e)  :\t", tf/opt->substeps);
    tic(0);

    /* Compute solution */
#if 1
    foreach_component(nc, P, NULL,
		      solveScalar,
		      solveSystem,
		      &sys, &rep);
#else
    solveSystem(0, n_cells, &sys, &rep);
#endif
    /* Swap u and u_old */
    double *tmp = sys.u; sys.u = sys.u_old; sys.u_old = tmp;

    /* Advance time */
    t  = t + sys.dt;
    sys.dt = sys.dt < tf - t ? sys.dt : tf - t;

    double elapsed = toc(0);
    Message(&str, &len, "Time elapsed: %10.7f seconds\n", elapsed);
    r->internal_time += elapsed;
  }

  /* ##### CHECK THIS! ######### */
  if (u == sys.u)
  {
    memcpy(u, sys.u_old, sys.n*sizeof(double));
  }

  /* Permute final solution back to original sequence */
  inverse_permute_dvector(u, work, rowP, n_cells);

  int i;
  int max=0;
  int sum=0;
  int num0=0;
  for (i=0; i<n_cells; ++i)
  {
    sum += rep.iterations[i];
    max  = max < rep.iterations[i] ? rep.iterations[i] : max;
    if (rep.iterations[i]==0) num0++;
  }

  Message (&str, &len, "Average iterations    :        %2.2f\n", sum*1.0/n_cells);
  Message (&str, &len, "Max iterations        :    %8d\n", max);
  Message (&str, &len, "Number of zero it.    :    %8d\n", num0);
  r->average_iterations = sum*1.0/n_cells;
  r->max_iterations     = max;
  r->num_zeroiterations = num0;
  loop_statistics(P, nc, &r->num_looped_cells, &r->num_loops, &r->max_loopsize);

  /* Free report */
  if (1)
  {
    free(rep.iterations);
    free(rep.type);
  }

  sparse_free(sys.V);
  free (rowP);
  free (P);
  free (work);

  return 0;
}







/*---------------------------------------------------------------------------*/
static int
insaneInput(const int nrhs, const mxArray *prhs[])
/*---------------------------------------------------------------------------*/
{
  if (nrhs != 6)
  {
    mexErrMsgTxt("Insane! Please call me with 6 input arguments\n");
    return 1;
  }

  int m = mxGetM(prhs[0]);
  int n = mxGetN(prhs[0]);

  /* Flux matrix */
  if (!mxIsSparse(prhs[0]) || (m!=n))
  {
    mexErrMsgTxt("First argument must be a square sparse matrix.\n");
    return 1;
  }
  /* Pore volume */
  if ((int)mxGetNumberOfElements(prhs[1]) != m || !isFull(prhs[1]))
  {
    mexErrMsgTxt(
	  "\nInput error:\n"
	  "Argument 2 (porevolume) must have same number of elements as\n"
	  "the number of rows in argument 1.\n"
	  );
    return 1;
  }
  /* Initial saturation */
  if ((int)mxGetNumberOfElements(prhs[2]) != m || !isFull(prhs[2]))
  {
    mexErrMsgTxt(
	  "\nInput error:\n"
	  "Argument 3 (initial saturation) must have same number of elements\n"
	  "as the number of rows in argument 1.\n"
	  );
    return 1;
  }
  /* Final time */
  if (mxGetNumberOfElements(prhs[3]) != 1 || !mxIsDouble(prhs[3]))
  {
    mexErrMsgTxt(
	  "\nInput error:\n"
	  "Argument 4 (final time) must be a double.\n"
	  );
    return 1;
  }
 /* Source term */
  if ((int)mxGetNumberOfElements(prhs[4]) != m || !isFull(prhs[4]))
  {
    mexErrMsgTxt(
	  "\nInput error:\n"
	  "Argument 5 (source term) must have same number of elements as\n"
	  "the number of rows in argument 1.\n"
	  );
    return 1;
  }
 /* Options */
  if ((int)mxGetNumberOfElements(prhs[5]) != 1 || !mxIsStruct(prhs[5]))
  {
    mexErrMsgTxt(
	  "\nInput error:\n"
	  "Argument 6 (options) must have one element and be a struct.\n"
	  );
    return 1;
  }

  return 0;
}


/*---------------------------------------------------------------------------*/
int isFull(const mxArray *a)
/*---------------------------------------------------------------------------*/
{
  return !(mxIsSparse(a) || mxIsStruct(a) || mxIsCell(a));
}

#if CHECK_DIVERGENCE
/*---------------------------------------------------------------------------*/
static int
checkDivergence(const sparse_t *A, const double *Q)
/*---------------------------------------------------------------------------*/
{


  int     m = A->m;
  int     i,j,k;
  double  val;

  double  *rowsum = calloc(m, sizeof(double));

  for(i=0; i<m; ++i)
  {
    for (k=A->ia [i]; k<A->ia [i+1]; ++k)
    {
      j   = A->ja [k];
      val = A->a  [k];

      if (val < 0.0)
      {
	mexErrMsgTxt ("Error!\n\n\tPositive fluxes expected.\n\n");

	free (rowsum);
	return 1;
      }
      rowsum [j] -= val;
      rowsum [i] += val;
    }
    rowsum [i] += Q[i];
  }

  double maxflux = 0.0;
  for (k=0; k<A->ia[m]; ++k)
  {
    if (fabs (A->a [k]) > maxflux)
    {
      maxflux = fabs (A->a [k]);
    }
  }


  double infnorm = 0.0;
  for(i=0; i<m; ++i)
  {
    if (fabs (rowsum [i]) > infnorm)
    {
      infnorm = fabs (rowsum [i]);
    }
  }

  free(rowsum);
  if (infnorm/maxflux > 1e-6)
  {
    mexErrMsgTxt("Error!\n\n\t||div(v)||/||v|| > 1e-6 (%e)\n\n", infnorm);
    return 1;
  }
  else
  {
    message(2, "OK\n");
  }

  return 0;
}
#endif
