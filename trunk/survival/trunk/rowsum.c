/*
**  SCCS $Id: rowsum.c,v 4.1 1993-06-29 10:06:53 therneau Exp $
**
** Add up data along rows
**
** Input
**      dim:   integer vector, the #rows and #columns of the matrix
**      x  :   matrix of data (remember, S uses column major order!)
**      group: the group to which each row belongs
**
** Output:
**      dd[0]: the number of unique groups found
**      x    : rows 1 to dd[0] contain the sums.
*/

void rowsum(dim, x, group)
long    *dim;
double  *x,
	*group;
    {
    register int i,j, k;
    int     nrow,
	    ncol;
    int     newrow;
    double  tgrp,
	    sum;
    double  dummy;

    nrow = dim[0];
    ncol = dim[1];

    dummy =0;
    for (i=0; i<nrow; i++) if (group[i] < dummy) dummy = group[i];
    dummy = (dummy/2) -1;    /*no group uses this number */

    newrow =0;
    for (i=0; i<nrow; i++) {
	if (group[i] > dummy) {
	    tgrp = group[i];
	    for (j=0; j<ncol; j++) {
		sum =0;
		for (k=i; k<nrow; k++)
		    if (group[k] == tgrp) sum += x[k + j*nrow];
		x[newrow + j*nrow] = sum;
		}
	    for (k=i; k<nrow; k++)
		if (group[k] == tgrp) group[k] = dummy;
	    newrow++;
	    }
	}
    dim[0] = newrow;
    }
