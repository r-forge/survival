/* SCCS $Id: pyears2.c,v 4.1 1993-12-02 21:36:43 therneau Exp $  */
/*
**  Person-years calculations.
**     same as pyears1, but no expected rates
**
**  Input:
**      n       number of subjects
**      ny      number of columns of y.
**      doevent does y have an 'events' column?  1=yes, 0=no
**              if ny=2 and doevent=1, then "start" is missing.
**      y[3,n]  contains start, stop, and event for each subject
**
**  output table's description
**      odim        number of dimensions
**      ofac[odim]  1=is a factor, 0=continuous (time based)
**      odims[odim] the number of rows, columns, etc
**      ocut[]      for each non-factor dimension, the odim[i]+1 cutpoints
**                        that define the intervals; concatonated.
**      odata[odim, n]  the subject data-- where each indexes into the
**                        expected table, at time 0.
**
** Output:
**      pyears     output table of person years
**      pn         number of observations that contribute to each cell
**      pcount     number of events
**      offtable   total person years that did not fall into the output table
**
** Scratch     allocated on the fly
**      scratch[edim]
*/

double **dmatrix();
double pystep();

/* names that begin with "s" will be re-declared in the main body */
void pyears2(sn, sny, sdoevent, sy,
		     sodim, ofac, odims, socut, sodata,
		     pyears, pn, pcount, offtable)

long    *sn,
	*sny,
	*sdoevent,
	*sodim,
	ofac[],
	odims[];

double  *sy,
	*socut,
	*sodata,
	*pyears,
	*pn,
	*pcount,
	*offtable;
    {
    register int i,j;
    int     n,
	    ny,
	    doevent,
	    odim;
    double  *start,
	    *stop,
	    *event,
	    **ocut,
	    **odata;
    double  *data;
    double  timeleft,
	    thiscell;
    int     index;
    int     dostart;
    int     d1;    /* a dummy */
    double  d2;    /* a dummy for pystep */

    n = *sn;
    ny= *sny;
    doevent = *sdoevent;
    odim = *sodim;
    start = sy;
    if (ny==3 || (ny==2 && doevent==0)) {
	stop = sy +n;
	dostart =1;
	}
    else   {
	stop  = sy;
	dostart =0;
	}
    event = stop +n;
    odata = dmatrix(sodata, n, odim);
    data  = (double *) S_alloc(odim, sizeof(double));
    /*
    ** will be a ragged array
    */
    ocut = (double **)S_alloc(odim, sizeof(double *));
    for (i=0; i<odim; i++) {
	ocut[i] = socut;
	if (ofac[i]==0) socut += odims[i] +1;
	}

    *offtable =0;
    for (i=0; i<n; i++) {
	/*
	** initialize
	*/
	for (j=0; j<odim; j++) {
	    if (ofac[j] ==1 || dostart==0) data[j] = odata[j][i];
	    else                           data[j] = odata[j][i] + start[i];
	    }
	if (dostart==1) timeleft = stop[i] - start[i];
	else            timeleft = stop[i];
	/*
	** add up p-yrs
	*/
	while (timeleft >0) {
	    thiscell = pystep(odim, &index, &d1, &d2, data, ofac, odims, ocut,
				    timeleft, 0);
	    if (index >=0) {
		pyears[index] += thiscell;
		pn[index] += 1;
		}
	    else *offtable += thiscell;

	    for (j=0; j<odim; j++)
		if (ofac[j] ==0) data[j] += thiscell;
	    timeleft -=thiscell;
	    }
	if (index >=0 && doevent) pcount[index] += event[i];
	}
    }
