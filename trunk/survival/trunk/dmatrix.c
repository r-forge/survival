/* SCCS $Id: dmatrix.c,v 4.4 1992-08-25 14:32:29 grill Exp $  */
/*
** set up ragged arrays, with #of columns and #of rows
*/
double **dmatrix(array, ncol, nrow)
double  *array;
int ncol, nrow;
    {
    register int i;
    register double **pointer;

    pointer = (double **) S_alloc(nrow, sizeof(double *));
    for (i=0; i<nrow; i++) {
	pointer[i] = array;
	array += ncol;
	}
    return(pointer);
    }
