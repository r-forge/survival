/*
**  SCCS @(#)survdiff.c	2.1 6/12/91
*/
#include <math.h>

survdiff2(nn, nngroup, rho, time, status, group, obs, exp, var, risk, kaplan)
long    *nn,
	*nngroup,
	group[],
	status[];
double  *rho;
double  time[],
	obs[],
	exp[],
	var[],
	risk[],
	kaplan[];
    {
    register int i,j,k;
    int kk;
    int n, ngroup;
    double km, nrisk, wt, tmp;
    double deaths;

    n = *nn;
    ngroup = *nngroup;

    /*
    ** Compute the k-m, which is only needed if rho!=0
    **   We want it set up as a left-continuous function (unusual)
    */
    if (*rho !=0){
	km =1;
	for (i=0; i<n; ) {
	    nrisk = n-i;
	    deaths =status[i];
	    for (j=i+1; time[j]==time[i] & j<n; j++) {
		deaths += status[j];
		}
	    kaplan[j-1] = km;
	    km = km * (nrisk-deaths)/nrisk;
	    i=j;
	    }
	}

    /*
    ** Now for the actual test
    */
    for (i=0; i<ngroup; i++) {
	obs[i]=0;
	exp[i]=0;
	risk[i]=0;
	}
    for (i=0; i< (ngroup-1)*(ngroup-1); i++)  var[i]=0;

    for (i=n-1; i>=0; i--) {
	deaths = 0;
	if (*rho ==0) wt=1;
	else          wt= pow(kaplan[i], *rho);
	for (j=i; time[j]==time[i] & j>=0; j--) {
	    k = group[j]-1;
	    deaths += status[j];
	    risk[k] += 1;
	    obs[k] += status[j] *wt;
	    }
	i=j +1;
	nrisk = n-i;

	for (j=0; j<ngroup; j++) exp[j] += wt* deaths * risk[j] / nrisk;

	if (nrisk==1) continue;  /*only 1 subject, so no variance */
	kk =0;
	wt = wt*wt;
	for (j=1; j<ngroup; j++) {
		tmp = wt* deaths* risk[j]* (nrisk-deaths)/(nrisk *(nrisk-1));
		var[kk+j-1] += tmp;
		for (k=1; k<ngroup; k++) {
		    var[kk] -= tmp * risk[k] / nrisk;
		    kk++ ;
		    }
		}
	}
    }
