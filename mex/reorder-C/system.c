#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#include "sparse.h"
#include "nlsolvers.h"
#include "system.h"




#if MEXFILE
#include <mex.h>
extern int interrupted;
#define print   mexPrintf
#define malloc  mxMalloc
#define calloc  mxCalloc
#define realloc mxRealloc
#define free    mxFree
#endif



static void compute_upwind_flux(sparse_t*, double*, double*, sparse_t**, double *work);

static void compute_residual   (system_t *sys, int begin, int end,
			 double *u, double *residual);

static void compute_jacobian   (system_t *sys, int begin, int end,
			 double *u, sparse_t *jacobian, double *work);

/* Error_1 Not block lower triangular */
static const char err_100 [] =
  "system::compute_jacobian assumes sys->V to be block lower triangular\n"
  "with respect to arguments begin and end.  Here, begin = %d, \n"
  "end = %d and column %d was found.";



/* visc_w/visc_o */
double mobratio = 1.0/1.0;
/*---------------------------------------------------------------------------*/
static double fluxfun(double sw)
/*---------------------------------------------------------------------------*/
{
  sw = sw > 1.0 ? 1.0 : sw;
  sw = sw < 0.0 ? 0.0 : sw;

  double so = 1.0-sw;
  return sw*sw/(sw*sw + mobratio*so*so);
}

/*---------------------------------------------------------------------------*/
static double dfluxfun(double sw)
/*---------------------------------------------------------------------------*/
{
  sw = sw > 1.0 ? 1.0 : sw;
  sw = sw < 0.0 ? 0.0 : sw;

  double so    =  1.0-sw;
  double denom = (sw*sw + mobratio*so*so);
  denom *= denom;
  return 2.0*mobratio*so*sw/denom;
}





/*
 * Input arguments:
 *
 *     vec      - pointer to double vector
 *
 *     n        - size of vec
 *
 * Output:
 *
 *     return   - inf norm of vec
 *
 */
/*---------------------------------------------------------------------------*/
static double
infnorm(double *vec, const int n)
/*---------------------------------------------------------------------------*/
{
  int    j;
  double norm = 0.0;
  for (j=0; j<n; ++j)
  {
    double val = fabs(vec[j]);
    if (val > norm)norm = val;
  }
  return norm;
}



static const char warning_str[]=
  "  In cell %d, the target function G(s) = s+Af(s)-B cannot have a zero in\n"
  "  [0,1] if either A + 1 < B or B < 0.  In this case A=%5f and B=%5f.\n";


static void   extractScalarSubProblem(int, system_t*, double*);


/*
 *
 *
 *
 * Without gravity, f is strictly increasing with range [0,1].
 * Consequently,
 *
 *   G = S + Af(s) -B
 *
 * is also strictly increasing.  Furthermore,
 *
 *   G(0) = -B
 *
 *   G(1) = 1 + A - B
 *
 * Therefore, G=0 has a solution in [0,1] if  1+A > B > 0.
 */




/*---------------------------------------------------------------------------*/
static double G(double s, void *data)/* double A, double B) */
/*---------------------------------------------------------------------------*/
{
  double *D = (double *) data;
  return s + D[0]*fluxfun(s) - D[1];
}

/*---------------------------------------------------------------------------*/
static double dG (double s, void *data)/* double A) */
/*---------------------------------------------------------------------------*/
{
  double *D = (double *) data;
  return 1.0 + D[0]*dfluxfun(s);
}


/*Question:  What about residual oil or water saturations?  */
/*
 * d   - data    (flux, pore volume, old saturations, ...)
 * r   - report  (number of iterations, ...)
 * opt - options (tolerance, max iteratuions, solver type, ...)
 */
/*---------------------------------------------------------------------------*/
int solveScalar(int begin, int end, void *d, void *r)
/*---------------------------------------------------------------------------*/
{
  system_t *sys = d;
  report_t *rep = r;
  assert(end-begin == 1);
  int cell = begin;

  enum Method {NEWTON, BRENT, RIDDER, REGULAFALSI, BISECTION} method = REGULAFALSI;
  double       TOL    = sys->tolerance;/* 1.0e-6; */
  const int    maxit  = sys->maxit;
  double      *s      = sys->u;

  double       data[2];
  extractScalarSubProblem(cell, sys, data);

  /* Check if solution in [0,1] exist */
  if (( data[0] + 1 < data[1] ) ||  ( data[1] < 0 ))
  {

#if 1
    if    ( data[1] <= 0.0)         s [cell] = 0.0;
    else /* data[1] >= data[0]+1 */ s [cell] = 1.0;
#else
    print(warning_str, cell, data[0], data[1]);
#endif

  }
  else
  {
    switch (method)
    {
      default:
      case BISECTION:
	s [cell] = bisection  ( &G, data, TOL, maxit, &rep->iterations[cell]);
	break;

      case BRENT:
	s [cell] = brent      ( &G, data, TOL, maxit, &rep->iterations[cell]);
	break;

      case RIDDER:
	s [cell] = ridder     ( &G, data, TOL, maxit, &rep->iterations[cell]);
	break;

      case REGULAFALSI:
	s [cell] = regulafalsi( &G, data, TOL, maxit, &rep->iterations[cell]);
	break;

      case NEWTON:
	s [cell] = newton     ( &G, &dG, data, TOL, maxit, &rep->iterations[cell]);
	break;
    }
  }

  sys->f[cell]   = fluxfun(s[cell]);
  rep->type[cell] = 0;
  return 0;
}


/*---------------------------------------------------------------------------*/
static void
extractScalarSubProblem(int            cell,
			system_t      *sys,
			double        *D)
/*---------------------------------------------------------------------------*/
{
  int k;
  double sumIn  = 0.0;
  double sumOut = 0.0;

  for (k=sys->V->ia[cell]; k<sys->V->ia[cell+1]; ++k)
  {
    int j = sys->V->ja[k];
    assert ( j <= cell );

    if (j<cell)
    {
      /* Incoming flux */
      sumIn -= sys->V->a[k]*sys->f[j];
    }
    else
    {
      /* Outgoing flux and negative source */
      sumOut   += sys->V->a[k];
    }
  }

  /* Incoming source */
  sumIn -= sys->Qp[cell]*fluxfun(1.0);


  D[0]   = sumOut *sys->dt;
  D[1]   = sumIn  *sys->dt + sys->u_old[cell];
}




/*
 * Solve
 *
 *     s^n - s^{n-1} + dt*flux*f(s^n) = max(sources,0)*f(s^n)
 *
 * using an implicit upwind scheme.  Assume that flux has a block lower
 * triangular sparsity structure, such that s_i^n for i<cell is already
 * computed, and the flux function f(s^n) is correct.
 *
 * The nonlinear iteration use two different strategies to ensure that a
 * solution is found eventually.  If convergence is not achieved within
 * a certain number of iterations, the time step is halved.  For each
 * time step, line search is used to improve robusteness.
 *
 * Input arguments:
 *
 *     cell     -
 *
 *     sz       -
 *
 *     sat      - s^n
 *
 *     sato     - s^{n-1}
 *
 *     ff       - f(s^n)
 *
 *     flux     -
 *
 *     Q        - max(sources, 0)
 *
 *     DT       -
 *
 * Output:
 *
 *     sat      -
 *
 *     sato     -
 *
 *     ff       -
 *
 */

/*---------------------------------------------------------------------------*/
int solveSystem(int begin, int end, void *d, void *r)
/*---------------------------------------------------------------------------*/
{
  system_t *sys = d;
  report_t *rep = r;

/*   double t = 0.0; */
  int    i;
  double TOL   = sys->tolerance;/* 1.0e-8; */
  int    maxit = sys->maxit;

  int    size  = end-begin;
  int    nnz   = sys->V->ia[end]-sys->V->ia[begin];

  /* external fluxes and sources */
  double *G    = malloc (  size * sizeof (double));
  double *ds   = malloc (  size * sizeof (double));
  double *work = malloc (5*sys->V->m * (sizeof (double)));
  sparse_t *dG = sparse_alloc (size, size, nnz);



  double onorm = TOL+1.0;
  int    it    = 0;

  /* Neutral initial guess */
  for (i=begin; i<end; ++i)
  {
    sys->u[i] = 0.5;
  }

  /* Initial residual */
  compute_residual   (sys, begin, end, sys->u, G);
  onorm = infnorm(G, size);

  while ((onorm > TOL) && (it++ < maxit) && !interrupted)
  {
    compute_jacobian   (sys, begin, end, sys->u, dG, work);
    sparse_solve(dG, G, ds);

    for (i=begin; i<end; ++i)
    {
      sys->u [i] -= 0.8*ds[i-begin];
    }


   /* Clip solution */
    for (i=begin; i<end; ++i)
    {
      sys->u[i] = (sys->u[i] > 1.0? 1.0 : sys->u[i]);
      sys->u[i] = (sys->u[i] < 0.0? 0.0 : sys->u[i]);
    }
    compute_residual   (sys, begin, end, sys->u, G);
    onorm = infnorm(G, size);
  }

  for (i=begin; i<end; ++i)
  {
    rep->iterations[i] = it;
    rep->type[i] = 1;
  }

  free (work);
  free (G);
  free (ds);
  sparse_free(dG);
  return 0;
}















/*
 * - For better memory management, a sparse matrix of fixed size should be able to
 *   borrow storage from work array.
 *
 *
 *
 */

/*---------------------------------------------------------------------------*/
void init_system        (system_t *sys,
			 double    dt,
			 sparse_t *flux,
			 double   *sources,
			 double   *volumes,
			 double   *u,
			 double    mr,
			 int       maxit,
			 double    tolerance,
			 double   *work)
/*---------------------------------------------------------------------------*/
{
  /*  work must be 2n * sizeof(double) bytes long.  Work is here to avoid
   *  malloc/free inside this function.   */

  size_t     n          = flux->m;
  sparse_t *upwind_flux = {0};

  /* Allocates a new sparse matrix. */
  compute_upwind_flux(flux, sources, volumes, &upwind_flux, work);

  sys->n      = n;
  sys->dt     = dt;
  sys->V      = upwind_flux;

  sys->u      = u;
  sys->u_old  = work;
  sys->f      = work +   n;
  sys->Qp     = sources;

  memcpy(sys->u_old, sys->u, sys->n*sizeof(double));


  sys->maxit = maxit;
  sys->mr = mr;
  sys->tolerance = tolerance;

#if 0

  sys->A      = A;
  sys->B      = B;

  /* Scale pressure terms */
  if (sys->A)
  {
    size_t i;
    for (i=0; i<n; ++i)
    {
      sys->A[i] /= volumes[i];
    }
  }

  if (sys->B)
  {
    int i;
    for (i=0; i<n; ++i)
    {
      sys->B[i] /= volumes[i];
    }
  }
#endif
}


/* (1+dtA)u^n + dt(B+V)f(u^n) = dtQ +u^{n-1} */
/*---------------------------------------------------------------------------*/
static void compute_residual   (system_t *sys, int begin, int end, double *u, double *residual)
/*---------------------------------------------------------------------------*/
{
  int i,k;

  /* We assume that sys->f has been filled with
   * fluxfun(u[0:begin-1]) already    */
  for (i=begin; i<end; ++i)
  {
    sys->f[i]  =  fluxfun(u[i]);
  }

  /* Clear end of sys->f.  We could assume f = 0 for
   * all i > end, but this would be a bug waiting to bite */
  for(i=end; i<sys->n; ++i)
  {
    sys->f[i] = 0.0;
  }

  for (i=begin; i<end; ++i)
  {
    /* Compute terms (1 + dt*A) *u - u_old + dt*Qp */
    residual [ i-begin ] = (/*(1.0 + sys->dt * sys->A [i]) * */
			    u [i] - sys->u_old [i]   +   sys->dt * sys->Qp [i]);


    /*  Compute term dt*(V+B)*f(u) */
    int b = sys->V-> ia [i];
    int e = sys->V-> ia [i+1];
    for (k=b; k<e; ++k)
    {
      const double v = sys -> V->a [k];/* + sys->B[i];*/
      const int    j = sys -> V->ja[k];
      residual [ i-begin ] += sys->dt * v * sys->f[ j ];
    }
  }
}

/* work must have room for (end-begin)* sizeof double bytes. */
/*---------------------------------------------------------------------------*/
static void compute_jacobian   (system_t *sys, int begin, int end, double *u, sparse_t *jacobian, double *work)
/*---------------------------------------------------------------------------*/
{
  assert (begin >= 0     );
  assert (end   <= sys->n);
  int i,k;
  double *df = work;
  /* Assume jacobian has the same structure as the block V(begin:end-1, begin:end-1)*/

  /* Clear jacobian */
  sparse_zeros (jacobian);

  /* Compute derivative of flux function in grid blocks */
  for (i=begin; i<end; ++i)
  {
    df[i] = dfluxfun (u [i]);
  }

  int pos = 0;
  jacobian->ia[0] = 0;
  for (i=begin; i<end; ++i)
  {

    int b = sys->V-> ia [i];
    int e = sys->V-> ia [i+1];
    for (k=b; k<e; ++k)
    {
      assert (pos <= jacobian->nres);
      int j = sys->V->ja[k];

      if (j < begin)
      {
	; /* do nothing*/
      }
      else if (j < end)
      {
	jacobian->ja [pos] = j-begin;
	jacobian->a  [pos] = sys->dt * ( sys->V->a[k]/*  + sys->B[i] */ ) * df[i]; /* ??? */
	if (i==j)
	{
	  jacobian->a[pos] += 1.0;/*  + sys->dt * sys->A[i]; */
	}
	++pos;
      }
      else
      {
	print(err_100, begin, end, j);
      }
      jacobian->ia[i-begin+1]=pos;
    }
  }
}


/*---------------------------------------------------------------------------*/
static void
compute_upwind_flux(sparse_t *Fin, double *source, double *volume,
		    sparse_t **upwind_flux, double *work)
/*---------------------------------------------------------------------------*/
{
  /* work must be n*sizeof(double) bytes long. */
  int       i,k;
  const int m = Fin->m;
  const int N = Fin->ia[m];

  sparse_t *F = sparse_alloc(m, m, N+m);

  /* Copy Fin to F and divide by PV */
  for (i=0; i<m+1; ++i)
  {
    F->ia[i] = Fin->ia[i]+i;
  }

  for (i=0; i<m; ++i)
  {
    for (k=Fin->ia[i]; k<Fin->ia[i+1]; ++k)
    {
      F->ja[k+i] =  Fin->ja[k];
      F->a [k+i] = -Fin->a [k]/volume[i];
    }
  }


  for (i=0; i<N; ++i)
  {
    work[i] = 0.0;
  }


  /* Add outgoing fluxes */
  double *r = Fin->a;
  int    *q = Fin->ja;
  for (i=0; i<N; ++i)
  {
    work[*q++] += *r++;
  }

  /* Add negative source */
  for (i=0; i<m; ++i)
  {
    work  [i] += source[i] < 0 ? -source[i]               : 0.0;
    source[i]  = source[i] > 0 ? -source[i]*1.0/volume[i] : 0.0;
  }


  /* Put outgoing fluxes on diagonal. */
  for (i=0; i<m; ++i)
  {
    k = F->ia[i+1]-1;

    F->ja[k] = i;
    F->a [k] = work[i]/volume[i];

  }

  *upwind_flux = F;
}


#if 0

/* Solve whole system */
/*---------------------------------------------------------------------------*/
void
implicit(system_t *sys)
/*---------------------------------------------------------------------------*/
{
  double t = 0.0;
  int i,j;
  int          it    = 0;
  double       relax = 1.0;
  double       TOL   = 1.0e-3;
  int          maxit = 2500;

  double dt          = sys->dt;

  ssystem_t S = {sys->A, sys->Q, sys->Sold, sys->S};
  double  *G;
  double  *ds;
  double  *df;


  df      = (double *)malloc ( size *     sizeof (double ));
  G       = (double *)malloc ( size *     sizeof (double)  );
  ds      = (double *)malloc ( size *     sizeof (double)  );

  message(3, "\nTolerance          :\t%e\n", TOL);
  message(3, "------------------------------------------------------\n");
  message(3, " Time interval           it  residual     relax       \n");
  message(3, "------------------------------------------------------\n");
  tic(4);

  /*   extractSparseSubProblem(begin, end, sys, &S); */
  assert(S.A->m==size);

  sparse_t *dG = NULL;
  copy_sparse (&dG, S.A);

   /* memcpy(S.xold, sys->Sold, size*sizeof(double)); */


  double onorm;
  int    restart;
  while ((t < sys->dt) && (!interrupted))
  {

    restart = 1;
    while(restart && !interrupted)
    {
      onorm = TOL+1.0;
      it    = 0;
#if 1
      /* Neutral initial guess */
      for (j=0; j<size; ++j) S.x[j] = 0.5;
#else
      /* Neutral initial guess */
      for (j=0; j<size; ++j) S.x[j] = S.xold[j];
#endif
      /* Initial residual */
      compValue (&S, dt, sys->f, G);
      while ((onorm > TOL) && (it++ < maxit) && !interrupted)
      {
        compJacobian(&S, dt, df, dG);

        sparse_solve_block_triangular(dG, G, ds, ncomp, comp);
/*      sparse_solve(dG, G, ds); */

        relax = 0.9;

        doLinesearch(&S, ds, dt, &relax, &onorm);
        memcpy(G, ds, size*sizeof(double));

        /* Clip solution */
        for (i=0; i<size; ++i)
        {
          S.x[i] = (S.x[i] > 1.0? 1.0 : S.x[i]);
          S.x[i] = (S.x[i] < 0.0? 0.0 : S.x[i]);
        }
        message(3, "[%f %f]  %d   %f   %f  ", t, t+dt, it, onorm, relax);
        toc(4);tic(4);
      }

      if (onorm > TOL)  dt      = dt/2.0;
      else              restart = 0;

    }
    t  = t + dt;
    dt = (dt > sys->dt-t ? sys->dt-t : dt);
  }

#if 0
  /* Check output */
  for (i=0; i<size; ++i)
  {
    if (S.x[i] > 1.0) print("Warning! s[%d] > 1.0\n", i);
    if (S.x[i] < 0.0) print("Warning! s[%d] < 0.0\n", i);

  }
#endif

  free (df);
  free (G);
  free (ds);
  free_sparse(dG);
}
#endif
/*---------------------------------------------------------------------------*/
int implicit(void *d, void *r)
/*---------------------------------------------------------------------------*/
{
  system_t *sys = d;
  report_t *rep = r;

  double t = 0.0;
  int    i;
  double TOL   = sys->tolerance;/* 1.0e-8; */
  int    maxit = sys->maxit;

  int begin =0;
  int end = sys->V->m;
  int    size  = end-begin;
  int    nnz   = sys->V->ia[end]-sys->V->ia[begin];

  /* external fluxes and sources */
  double *G    = malloc (  size * sizeof (double));
  double *ds   = malloc (  size * sizeof (double));
  double *work = malloc (5*sys->V->m * (sizeof (double)));
  sparse_t *dG = sparse_alloc (size, size, nnz);
  double DT = sys->dt;

  int    it;
  int    restart;
  while ((t < DT) && (!interrupted))
  {

    restart = 1;
    while(restart && !interrupted)
    {

      double onorm = TOL+1.0;
      it    = 0;

      /* Neutral initial guess */
      for (i=begin; i<end; ++i)
      {
	sys->u[i] = 0.5;
      }

      /* Initial residual */
      compute_residual   (sys, begin, end, sys->u, G);
      onorm = infnorm(G, size);

      while ((onorm > TOL) && (it++ < maxit) && !interrupted)
      {
	compute_jacobian   (sys, begin, end, sys->u, dG, work);
	sparse_solve(dG, G, ds);

	for (i=begin; i<end; ++i)
	{
	  sys->u [i] -= 0.8*ds[i-begin];
	}


	/* Clip solution */
	for (i=begin; i<end; ++i)
	{
	  sys->u[i] = (sys->u[i] > 1.0? 1.0 : sys->u[i]);
	  sys->u[i] = (sys->u[i] < 0.0? 0.0 : sys->u[i]);
	}
	compute_residual   (sys, begin, end, sys->u, G);
	onorm = infnorm(G, size);
      }
      if (onorm > TOL)  sys->dt      = sys->dt/2.0;
      else              restart = 0;

    }
    t  = t + sys->dt;
    sys->dt = (sys->dt > DT-t ? DT-t : sys->dt);
  }

  for (i=begin; i<end; ++i)
  {
    rep->iterations[i] = it;
    rep->type[i] = 1;
  }
  sys->dt = DT;
  free (work);
  free (G);
  free (ds);
  sparse_free(dG);
  return 0;
}
