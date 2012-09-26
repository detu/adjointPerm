/*------------------------------------------------------------
 File       : strongcomp.h
 Created    : Wed Aug  6 10:42:05 2008
 Author     : Jostein R. Natvig <Jostein.R.Natvig@sintef.no>
 Revision   : $Id$


 Description
 Tarjans algorithm for finding strong components and topological
 sequence of a directed graph.  The graph is encoded in compressed
 sparse row format.

 Changes

------------------------------------------------------------*/

#ifndef STRONGCOMP_H
#define STRONGCOMP_H
void strongcomp(int size, int *ia, int *ja, int *rowP, int *P, int *ncomp, int *work);
int reorder   (sparse_t *A, int* rowP, int *P, void *work);
int foreach_component(int   nc,                    /* num components            */
	              int  *comp,                  /* begin and end of comp     */
		      int  *order,                 /* sequence                  */
		      int (*solvescalar)(int, int, void *, void*), /* one cell  */
		      int (*solvesystem)(int, int, void *, void*), /* may cells */
		      void *data,
		      void *report);
#endif /* STRONGCOMP_H */
