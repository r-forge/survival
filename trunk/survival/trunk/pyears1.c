/* SCCS $Id: pyears1.c,v 4.1 1993-12-02 21:36:42 therneau Exp $  */
/*
**  Person-years calculations, in its most general
**
**  Input:
**      n       number of subjects
**      ny      number of columns of y.
**      doevent does y have an 'events' column?  1=yes, 0=no
**              if ny=2 and doevent=1, then "start" is missing.
**      y[3,n]  contains start, stop, and event for each subject
**
**    expected table
**      edim        number of dimensions of the expected table
**      efac[edim]  1=is a factor, 0=continuous (time based)
**                    >=2  special handling for US "calendar year"
**      edims[edim] the number of rows, columns, etc
**      ecut[ ]     the starting points for each non-factor dimension,
**                          strung together.
**      expect      the actual table of expected rates
**      edata[edim, n]  the subject data-- where each indexes into the
**                        expected table, at time 0.
**
**   output table's description
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
**      pexpect    expected number of events
**      offtable   total person years that did not fall into the output table
**
** Scratch   -- allocated on the fly
**      scratch[edim + odim]
*/

double **dmatrix();
double pystep();

/* names that begin with "s" will be re-declared in the main body */
void pyears1(sn, sny, sdoevent, sy,
		     sedim, efac, edims, secut, expect, sedata,
		     sodim, ofac, odims, socut, sodata,
		     pyears, pn, pcount, pexpect, offtable)

long    *sn,
	*sny,
	*sdoevent,
	*sedim,
	*sodim,
	efac[],
	ofac[],
	edims[],
	odims[];

double  *sy,
	*secut,
	*socut,
	*expect,
	*sedata,
	*sodata,
	*pyears,
	*pn,
	*pcount,
	*pexpect,
	*offtable;
    {
    register int i,j;
    int     n,
	    ny,
	    doevent,
	    edim,
	    odim;
    double  *start,
	    *stop,
	    *event,
	    **ecut,
	    **ocut,
	    **edata,
	    **odata;
    double  *data,
	    *data2;
    double  timeleft,
	    thiscell,
	    etime,
	    et2;
    int     index,
	    indx, indx2;
    double  wt;
    int     dostart;

    n = *sn;
    ny= *sny;
    doevent = *sdoevent;
    edim = *sedim;
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
    edata = dmatrix(sedata, n, edim);
    odata = dmatrix(sodata, n, odim);
    i=edim + odim;
    data  = (double *) S_alloc(i, sizeof(double));
    data2 = data + odim;
    /*
    ** ecut and ocut will be ragged arrays
    */
    ecut = (double **)S_alloc(edim, sizeof(double *));
    for (i=0; i<edim; i++) {
	ecut[i] = secut;
	if (efac[i]==0)     secut += edims[i];
	else if(efac[i] >1) secut += (edims[i]*efac[i] * 1 - efac[i]);
	}

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
	for (j=0; j<edim; j++) {
	    if (efac[j] ==1 || dostart==0) data2[j] = edata[j][i];
	    else                           data2[j] = edata[j][i] + start[i];
	    }
	if (dostart==1) timeleft = stop[i] - start[i];
	else timeleft= stop[i];

	/*
	** add up p-yrs
	*/
	while (timeleft >0) {
	    thiscell = pystep(odim, &index, &indx2, &wt, data, ofac, odims,
				     ocut, timeleft, 0);
	    if (index >=0) {
		pyears[index] += thiscell;
		pn[index] += 1;

		/* expected calc */
		etime = thiscell;
		while (etime >0) {
		    et2 = pystep(edim, &indx, &indx2, &wt, data2, efac,
				 edims, ecut, etime, 1);
		    if (wt <1) pexpect[index]+= et2*(wt*expect[indx] +
						     (1-wt)*expect[indx2]);
		    else       pexpect[index]+= et2* expect[indx];
		    for (j=0; j<edim; j++)
			if (efac[j] !=1) data2[j] += et2;
		    etime -= et2;
		    }
		}
	    else  *offtable += thiscell;

	    for (j=0; j<odim; j++)
		if (ofac[j] ==0) data[j] += thiscell;
	    timeleft -=thiscell;
	    }
	if (index >=0 && doevent) pcount[index] += event[i];
	}
    }
