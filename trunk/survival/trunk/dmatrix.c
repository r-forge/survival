/* SCCS $Id: dmatrix.c,v 4.2 1992-08-10 13:48:21 grill Exp $  */
/*
** set up ragged arrays, with #of columns and #of rows
*/
double **dmatrix2(array, ncol, nrow)
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
