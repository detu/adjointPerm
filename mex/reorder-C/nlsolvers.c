/*------------------------------------------------------------
 File       : nlsolvers.c
 Created    : Thu Oct 16 09:38:02 2008
 Author     : Jostein R. Natvig <Jostein.R.Natvig@sintef.no>
 Revision   : $Id$


 Description


 Changes

------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "nlsolvers.h"


#if MEXFILE
#include <mex.h>
extern int interrupted;
#define print   mexPrintf
#define malloc  mxMalloc
#define calloc  mxCalloc
#define realloc mxRealloc
#define free    mxFree
#endif


static const char no_root_str[]=
  "  In %s:\n"
  "  With G(0) =% 5f, G(1) =% 5f, G(s) cannot have a zero in [0,1]!\n";

static const double EPS = 1e-14;

/*---------------------------------------------------------------------------*/
double
brent (double (*G)(double, void*), void *data, double tol, int maxit, int *iterations)
/*---------------------------------------------------------------------------*/
{
  double sa, Ga;
  double sb, Gb;
  double sc, Gc;
  double s,  Gs;
  int it = 0;

  sa = 0.0; Ga = G(sa, data);
  sb = 1.0; Gb = G(sb, data);

  if (Ga*Gb > 0)
  {
    print(no_root_str, "brent", Ga, Gb);
    return -1.0;
  }

  /* swap? */
  double absGa = fabs(Ga);
  double absGb = fabs(Gb);
  if (absGa < absGb)
  {
    double tmp;
    tmp = Ga;
    Ga  = Gb;
    Gb  = tmp;

    tmp = sa;
    sa  = sb;
    sb  = tmp;
  }


  sc = sa;
  Gc = Ga;

  int mflag = 1;
  double d=sc;
  while( (fabs(Gb)>EPS) && (fabs(sb-sa) > tol) && (it++<maxit))
  {

    if ( (fabs(Ga-Gc)>EPS) && (fabs(Gb-Gc)>EPS) )
    {
      /* Inverse quadratic interpolation */
#if 0
      double R = Gb/Gc;
      double S = Gb/Ga;
      double T = Ga/Gc;
      double P = S*(T*(R-T)*(sc-sb) - (1-R)*(sb-sa));
      double Q = (R-1)*(S-1)*(T-1);
      s = sb + P/Q;
#else
      s = sa*Gb*Gc/((Ga-Gb)*(Ga-Gc))
	+ sb*Ga*Gc/((Gb-Ga)*(Gb-Gc))
	+ sc*Ga*Gb/((Gc-Ga)*(Gc-Gb)) ;
#endif
    }
    else
    {
      /* Secant rule */
      s = sb - Gb*(sb-sa)/(Gb-Ga);
    }



    double absssb  = fabs(s-sb);
    double abssbsc = fabs(sb-sc);
    double absscd  = fabs(sc-d);

    int t1 = ((3*sa+sb)/4.0 -s)*(sb -s) > 0;
    int t2 =   mflag  && absssb >= abssbsc/2.0;
    int t3 = (!mflag) && absssb >= absscd /2.0;

    if ((t1 || t2 || t3))
    {
      s = 0.5*(sa+sb);
      mflag = 1;
    }
    else
    {
      mflag = 0;
    }

    Gs = G(s, data);

    d = sc;

    sc = sb;
    Gc = Gb;

    if (Ga*Gs <= 0)
    {
      sb = s;
      Gb = Gs;
    }
    else
    {
      sa = s;
      Ga = Gs;
    }

    absGa = fabs(Ga);
    absGb = fabs(Gb);
    if (absGa < absGb)
    {
      double tmp;
      tmp = Ga;
      Ga  = Gb;
      Gb  = tmp;

      tmp = sa;
      sa  = sb;
      sb  = tmp;
    }
  }

  *iterations = it;
  return sb;
}



/* Start with bracket [0,1] with G(0)*G(1)<0.  */
/*---------------------------------------------------------------------------*/
double
ridder (double (*G)(double, void*), void *data, double tol, int maxit, int *iterations)
/*---------------------------------------------------------------------------*/
{
  *iterations = 0;
  /* Initial guess is interval [0,1] of course */
  double s0=0.0;
  double G0=G(s0, data);
  if (fabs(G0)<EPS) return s0;

  double s1=1.0;
  double G1=G(s1, data);
  if (fabs(G1)<EPS) return s1;

  if (G0*G1 > 0)
  {
    print(no_root_str, "ridder", G0, G1);
    return -1.0;
  }

  if (G0>G1)
  {
    double swap;
    swap = G0;
    G0   = G1;
    G1   = swap;
  }

  double s3=0;
  double G3=10;

  int it = 0;
  while ( (fabs(G3) > tol)  &&  (it++ < maxit))
  {
    /* find zero crossing of line segment [(s0,G0), (s1,G1)] */
    double s2   = 0.5*(s0+s1);
    double G2   = G(s2, data);
    double root = sqrt(G2*G2-G0*G1);
    if (fabs(root)<EPS)
      return -1.0; /* Hmmm */

    double sgn  = G0>G1 ? 1.0 : -1.0;
    s3          = s2 + ( s2-s0 )*sgn*G2/root;
    G3          = G(s3, data);


    /* if     G2*G3<0  */
    if (G2*G3 <= 0.0)
    {
      if (G2 > G3)
      {
	s0 = s3;
	G0 = G3;
	s1 = s2;
	G1 = G2;
      }
      else
      {
	s0 = s2;
	G0 = G2;
	s1 = s3;
	G1 = G3;

      }

    }
    else if (G0*G3 <= 0.0)
    {
      s1 = s3;
      G1 = G3;
    }
    else if (G1*G3 <= 0.0)
    {
      s0 = s3;
      G0 = G3;
    }
    else
    {
      print(
	    "In ridder:\n"
	    "G0=%10.10f, G1=%10.10f, G3=%10.10f\n",
	    G0, G1, G3
	    );
      getchar();
    }
  }

  *iterations = it;
  return s3;
}



/* Start with bracket [0,1] with G(0)*G(1)<0.  Search by finding zero crossing
   sN of line segment [(sL,GL), (sR,GR)].  Set SL=sN if G(sN<0), sR=sN
   otherwise.*/
/*---------------------------------------------------------------------------*/
double
regulafalsi (double (*G)(double, void*), void *data, double tol, int maxit,
	     int *iterations)
/*---------------------------------------------------------------------------*/
{
  *iterations = 0;

  /* Initial guess is interval [0,1] of course */
  double s0=0.0;
  double G0=G(s0, data);
  if (fabs(G0)<EPS) return s0;

  double s1=1.0;
  double G1=G(s1, data);
  if (fabs(G1)<EPS) return s1;

  if (G0*G1 > 0)
  {
    print(no_root_str, "regulafalsei", G0, G1);
    return -1.0;
  }

  if (G0>G1)
  {
    double swap;
    swap = G0;
    G0   = G1;
    G1   = swap;
  }

  double sn =-1.0; /* "Undefined" value */
  double Gn = 10;

  int it = 0;
  while ( (fabs(Gn) > tol)  &&  (it++ < maxit))
  {
    /* find zero crossing of line segment [(s0,G0), (s1,G1)] */
    sn = s0 - (s1-s0)*G0/(G1-G0);
/*     sn = (s0*G1 - s1*G0)/(G1-G0); */
    Gn = G(sn, data);

#if 0
    /* Unmodified Regula-Falsi */
    /* maintain bracket with G1>G0 */
    if ( Gn>0 )
    {
      G1 = Gn;
      s1 = sn;
    }
    else
    {
      G0 = Gn;
      s0 = sn;
    }
#else
    if ((Gn>0.0)==(G0>0.0))
    {
      s0 = s1;
      G0 = G1;
    }


    /* Modified Regula-Falsi*/
    else
    {
      /* Illinois method */
      const double gamma_illinois = 0.5;

      /* Pegasus method */
      const double gamma_pegasus  = G1/(G1+Gn);

      G0 *= gamma_pegasus;
    }
    s1 = sn;
    G1 = Gn;
#endif

  }

  *iterations = it;
  return sn;
}

/* Start with bracket [0,1] with G(0)*G(1)<0.  Search by finding sN=0.5(sL+sR).
   Set SL=sN if G(sN<0), sR=sN otherwise.*/
/*---------------------------------------------------------------------------*/
double
bisection (double (*G)(double, void*), void *data, double tol, int maxit,
	   int *iterations)
/*---------------------------------------------------------------------------*/
{
  *iterations = 0;

  /* Initial guess is interval [0,1] of course */
  double s0=0.0;
  double G0=G(s0, data);
  if (fabs(G0)<EPS) return s0;

  double s1=1.0;
  double G1=G(s1, data);
  if (fabs(G1)<EPS) return s1;

  if (G0*G1 > 0.0)
  {
    print(no_root_str, "bisection", G0, G1);
    return -1.0;
  }

  if (G0>G1)
  {
    double swap;
    swap = G0;
    G0   = G1;
    G1   = swap;
  }

  double sn;
  double Gn=1;

  int it=0;
  while ( (fabs(Gn)>tol) && (it++ < maxit) )
  {

    sn = 0.5*(s0+s1);
    Gn = G(sn, data);

    if ( Gn>0 )
    {
      G1 = Gn;
      s1 = sn;
    }
    else
    {
      G0 = Gn;
      s0 = sn;
    }
  }
  *iterations = it;
  if (it >= maxit) print("Warning: convergence criterion not met\n");
  return 0.5*(s0+s1);
}


/*---------------------------------------------------------------------------*/
double
newton (double (*G)(double, void*), double (*dG)(double, void*), void *data,
	double tol, int maxit, int *iterations)
/*---------------------------------------------------------------------------*/
{
  *iterations = 0;

  /* Initial guess is interval [0,1] of course */
  double s0=0.0;
  double G0=G(s0, data);

  double s1=1.0;
  double G1=G(s1, data);

  if (G0*G1 > 0.0)
  {
    print(no_root_str, "newton", G0, G1);
    return -1.0;
  }

  const double relax    = 0.9;
  double       residual = 1.0e10;
  double       Gn, dGn;
  double       sn       = 0.5*(s0+s1);
  int          it       = 0;

  while ((fabs(residual) > tol) && (it++ < maxit))
  {
    dGn = dG(sn, data);
    Gn  =  G(sn, data);

    residual = -Gn/dGn;
    sn      += relax*residual;
  }

  *iterations = it;
  if (it >= maxit)print("Max num iterations reached\n");
  return sn;
}
