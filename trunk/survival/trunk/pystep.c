/* SCCS $Id: pystep.c,v 4.1 1993-12-02 21:36:44 therneau Exp $  */
/*
** Returns the amount of time that will be spent in the current "cell",
**  along with the index of the cell (treating a multi-way array as linear).
** This is a basic calculation in all of the person-years work.
**
** Input
**      nc:         number of categories
**      data[nc]    start points, for the data values
**      fac[nc]     1: category is a factor, 0: it is continuous
**                  >=2: special handling for "years" dim of US rate tables
**      dims[nc]    the extent of each category
**      cuts[nc,dims+1] ragged array, containing the start for each interval
**      step        the amount of time remaining for the subject.
**      edge        if =0, then the cuts contain +1 obs, and we are strict
**                      about out-of-range cells.  If it is a 1, then the
**                      table is assummed to extend infinitly at the edges.
**
** Output
**      *index   linear index into the array
**      if *index == -1, then the returned amount of time is "off table";
**  if one of the dimensions has fac >1 --
**      *index2   second index for linear interpolation
**      *wt       a number between 0 and 1, amount of wt for the first index
**                this will be 1 if none of the dims have fac >1
**
** Return value     amount of time in indexed cell.
*/

double pystep(nc, index, index2, wt, data, fac, dims, cuts, step, edge)
int     nc,
	edge,
	*index,
	*index2;
long    fac[],
	dims[];
double  data[],
	**cuts,
	*wt,
	step;
    {
    register int i,j;
    double maxtime;
    double shortfall;
    double temp;
    int k, kk, dtemp;


    kk=1;
    *index =0;  *index2=0;
    *wt =1;
    shortfall =0;
    maxtime = step;
    for (i=0; i<nc; i++) {
	if (fac[i]==1) *index += (data[i]-1) * kk;
	else {
	    if (fac[i]>0) dtemp = fac[i]*dims[i];
	    else          dtemp = dims[i];
	    for (j=1; j<dtemp; j++) if (data[i] < cuts[i][j]) break;
	    if (edge==0 || j<dtemp) {
		temp = cuts[i][j] - data[i];
		if (temp <= 0) shortfall = step;
		else if (temp <maxtime)  maxtime = temp;
		}

	    j--;
	    if (fac[i]>1 ) {  /* return to actual indices, & interpolate */
		*wt = 1.0 - (j%fac[i])/ (double)fac[i];
		j /= fac[i];
		*index2 = kk;
		}

	    temp = cuts[i][j] - data[i];   /*if positive, we're < first cut */
	    if (edge==0 && temp > shortfall) {
		if (temp < step)  shortfall = temp;
		else              shortfall = step;
		}
	    else *index += j*kk;
	    }
	kk *= dims[i];
	}

    *index2 += *index;
    if (shortfall ==0) return(maxtime);
    else  {
	*index = -1;
	return(shortfall);
	}
    }
