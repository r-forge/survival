/*  SCCS $Id: dmatrix.c,v 5.1 1998-08-30 14:52:48 therneau Exp $
/*
** set up ragged arrays, with #of columns and #of rows
*/
#include "survproto.h"
#include "survS.h"

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
