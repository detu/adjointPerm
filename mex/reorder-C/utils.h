#ifndef _UTILS_H_
#define _UTILS_H_

void permute_dvector(double *vec, double *work, const int *p, const int n);
void inverse_permute_dvector(double *vec, double *work, const int *p, const int n);

void     tic(int timernum);
double   toc(int timernum);
#if 0
void     setVerboseLevel(const int level);

int      message (const int verbose, const char *format, ...);
int      warning (const char *format, ...);
int      error   (const char *format, ...);
#endif
#endif
