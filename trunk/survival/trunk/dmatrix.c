/* SCCS $Id: dmatrix.c,v 4.1 1992-03-04 16:51:50 therneau Exp $  */
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
