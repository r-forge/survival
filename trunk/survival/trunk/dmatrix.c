/*  SCCS $Id: dmatrix.c,v 5.2 1998-10-27 17:36:25 therneau Exp $
/*
** set up ragged arrays, with #of columns and #of rows
*/
#include "survS.h"
#include "survproto.h"

double **dmatrix(double *array, int ncol, int nrow)
    {
S_EVALUATOR
    register int i;
    register double **pointer;

    pointer = (double **) ALLOC(nrow, sizeof(double *));
    for (i=0; i<nrow; i++) {
	pointer[i] = array;
	array += ncol;
	}
    return(pointer);
    }
