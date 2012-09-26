#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <string.h>
#include <assert.h>

#include <mex.h>
#include "utils.h"
extern int interrupted;

static        int     verbose_level = 1;
static const  int     n_timers = 8;
static struct timeval tstart[8];


/*---------------------------------------------------------------------------*/
void
permute_dvector(double *vec, double *work, const int *p, const int n)
/*---------------------------------------------------------------------------*/
{
  int i;
  memcpy(work, vec, n * sizeof (double));
  for (i=0; i<n; ++i) vec[i] = work[p[i]];
}

/*---------------------------------------------------------------------------*/
void
inverse_permute_dvector(double *vec, double *work, const int *p, const int n)
/*---------------------------------------------------------------------------*/
{
  int i;
  memcpy(work, vec, n * sizeof (double));
  for (i=0; i<n; ++i) vec [p [i]] = work [i];
}



/*---------------------------------------------------------------------------*/
void tic (int timernum)
/*---------------------------------------------------------------------------*/
{
  if (timernum > n_timers)
  {
    error ("Timer number %d outside range\n", timernum);
    exit(1);
  }
  gettimeofday(&tstart[timernum], NULL);
}

/*---------------------------------------------------------------------------*/
double toc (int timernum)
/*---------------------------------------------------------------------------*/
{
  if (timernum > n_timers)
  {
    error ("Timer number %d outside range\n", timernum);
    exit(1);
  }
  struct timeval tend;
  gettimeofday(&tend, NULL);

  double seconds =
    difftime(tend.tv_sec, tstart[timernum].tv_sec)
    +1.0e-6*(tend.tv_usec-tstart[timernum].tv_usec);

/*   message(2, "Time elapsed: %10.7f seconds", seconds); */
}

#if 0
/*---------------------------------------------------------------------------*/
void setVerboseLevel (const int level)
/*---------------------------------------------------------------------------*/
{
  const int max_level = 10;
  if ((level >= 0) && (level <= max_level))
    verbose_level = level;
  else
    error("Verbosity level must be between 0 (silent) and %d (noisy).",
	  max_level);
}


/*---------------------------------------------------------------------------*/
int message (const int verbose, const char *format, ...)
/*---------------------------------------------------------------------------*/
{
  int retval = 0;
  if (verbose <= verbose_level) /* Allow variable verbosity */
  {
    const int size = 1024;
    char *message = mxMalloc(size*sizeof(char));

    va_list ap;
    va_start(ap, format);
    retval = vsnprintf(message, size, format, ap);
    va_end(ap);

    mexPrintf("%s", message);
    mxFree(message);
  }
  return retval;
}


/*---------------------------------------------------------------------------*/
int warning (const char *format, ...)
/*---------------------------------------------------------------------------*/
{
  int retval = 0;
  if (verbose_level >= 1) /* Suppress warnings in silent mode. */
  {
    const int size = 1024;
    char *message = mxMalloc(size*sizeof(char));
    int n;

    va_list ap;
    va_start(ap, format);
    n  = snprintf  (message,   size,  "Warning:\n");
    n += vsnprintf (message+n, size-n, format,  ap);
    va_end(ap);

    mexWarnMsgTxt(message);
    mxFree(message);
  }
  return retval;
}

/*---------------------------------------------------------------------------*/
int error (const char *format, ...)
/*---------------------------------------------------------------------------*/
{
  int retval = 0;
  if (1) /* Always show error messages */
  {
    const int size = 1024;
    char *message = mxMalloc(size*sizeof(char));
    int n;

    va_list ap;
    va_start(ap, format);
    n  = snprintf(message, size, "Error:\n");
    n += vsnprintf (message+n, size, format,  ap);
    va_end(ap);

    mexErrMsgTxt(message);
    mxFree(message);
  }
  return retval;
}
#endif
