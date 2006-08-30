/*  SCCS  $Id: surv_callback.c,v 1.3 2006-08-30 20:53:46 m015733 Exp $
** callback routines for the survreg "other" distributions
**  This is based on appendix A of Chambers
*/
#include "survS.h"
#include "survproto.h"

static Sint nframe;
static Sint n;
static s_object *expr1;   /* Code to be executed */
static double  *dptr;
/*
** The first routine sets things up
*/
s_object *init_survcall(s_object *frame, s_object *n2, s_object *expr) {
S_EVALUATOR
    s_object *eta; 
   	
    nframe = INTEGER_VALUE(frame);/* the frame number of the calling routine */
    n      = INTEGER_VALUE(n2);   /* remember the length of the object */
    eta    = NEW_NUMERIC(n);      /* create a vector to hold the values */
    ASSIGN_IN_FRAME("zzeta", eta, nframe); /*make it visible in parent */
    dptr   = NUMERIC_POINTER(eta);  /* this points to the data region */
    expr1  = expr;	           /* the expression to evaluate */
    return(frame);    /* the book implies that something must be returned */
    }

/*
** This part is called by the survreg5 function, to get the distribution
*/
void surv_callback(double *z, double *dist) {
S_EVALUATOR
    s_object *value;
    double  *vd;
    int i;
 
    /* plug in the values */
    for (i=0; i<n; i++) dptr[i] = z[i];

    /* evaluate the expression */
    value = AS_NUMERIC( EVAL_IN_FRAME(expr1, nframe));
    
    /* grab off the results */
    if (LENGTH(value) != n*5) PROBLEM "Wrong length for return value" ERROR;
    vd = NUMERIC_POINTER(value);
    for (i=0; i< n*5; i++) dist[i] = vd[i];
    }

