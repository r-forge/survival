/*
**  SCCS $Id: cox_callback.c,v 1.2 2001-01-19 15:41:15 therneau Exp $
** callback routines for the coxph frailty interface
**  This version has been updated to use the "Green book" macros
*/
#include "survS.h"
#include "survproto.h" 

static long nframe;
static s_object *coef1;  /*sparse frailties */
static s_object *coef2;  /* coefs for non-sparse penalized terms */
static s_object *expr1;  /*the evaluation expression for sparse */
static s_object *expr2;  /* the expression for non-sparse */
static double *cdata1;   /* pointer to the data portion of coef1 */
static double *cdata2;
static long length1;
static long length2;

/*
** The first routine just saves away the location of the calling frame,
**   allocates an object for the 'coef1' vector, assigns the vector
**   to the orignal frame, and saves the expression to be computed.
*/
s_object *init_coxcall1(s_object *frame, s_object *nfrail, s_object *expr) {
    S_EVALUATOR

    nframe = INTEGER_VALUE(frame);/* the frame number of the calling routine */
    length1= INTEGER_VALUE(nfrail);
    coef1  = NEW_NUMERIC(length1);
    ASSIGN_IN_FRAME("coef1", coef1, nframe); 
    cdata1 = NUMERIC_POINTER(coef1);
    expr1  = expr;
    return(frame);   /* From the docs, it appears that I have to 
				 return something */
    }

s_object *init_coxcall2(s_object *frame, s_object *ncoef, s_object *expr) {
    S_EVALUATOR

    nframe = INTEGER_VALUE(frame);/* the frame number of the calling routine */
    length2= INTEGER_VALUE(ncoef);
    coef2  = NEW_NUMERIC(length2);
    ASSIGN_IN_FRAME("coef2", coef2, nframe); 
    cdata2 = NUMERIC_POINTER(coef2);
    expr2  = expr;
    return(frame);
    }

/*
** This part is called by the coxfit4 function, to get the penalty terms
*/
void cox_callback (int which, double *coef, double *first, 
	           double *second, double *penalty, long *flag) {
    S_EVALUATOR
	int i, j;
    s_object *coxlist;
    s_object **lptr;    /* pointer to the list of elements in coxlist */
    double *dptr;
    long   *logical;
    long  n;

    /* 
    ** Assign the frailty coef in the parent frame,
    **  and evaluate the expression.
    ** It is important that coef1 and coef2 be treated as "read only" objects
    **  in the parent frame, since we replace (again and again) the elements
    **  of the vectors without checking to see that they have not been 
    **  supplanted by a new copy.  In some sense, calling ASSIGN_IN_FRAME
    **  each time we entered this routine would be safer, except for the fact
    **  that in this version (6.0) reassigning the same object into a frame
    **  more than once crashes Splus.
    ** The expression darned well better return a list.
    */
    if (which==1) {
	for (i=0; i<length1; i++) cdata1[i]= coef[i];
	coxlist = EVAL_IN_FRAME(expr1, nframe);
	n = length1;
	}
    else {
	for (i=0; i<length2; i++) cdata2[i]= coef[i];
	coxlist = EVAL_IN_FRAME(expr2, nframe);
	n = length2;
        }
    if (!IS_LIST(coxlist)) PROBLEM
			       "The expression expr%d did not return a list!",
			          which ERROR;
    
    /*
    ** Grab data back off of the list
    **  The C-code here is fragile: it assumes that the elements of the
    **  list are coef, first, second, penalty, and flag IN THAT ORDER
    */
    j = LENGTH(coxlist);
    if (j !=5) PROBLEM
	  "The expression expr%d returned a list of %d elements, %d required",
		   which, j  ERROR;
    lptr = LIST_POINTER(coxlist);

    dptr = NUMERIC_POINTER(lptr[0]);  /* coef, possibly recentered */
    j =    LENGTH(lptr[0]);
    if (j != n) PROBLEM "Wrong length for coef, want %d got %d", n, j ERROR;
    for (i=0; i<j; i++) coef[i] = dptr[i];

    dptr = NUMERIC_POINTER(lptr[1]);  /* first derivative */
    j =    LENGTH(lptr[1]);
    if (j != n) PROBLEM "Wrong length for first, want %d got %d", n, j ERROR;
    for (i=0; i<j; i++) first[i] = dptr[i];

    dptr = NUMERIC_POINTER(lptr[2]);  /* second derivative */
    j =    LENGTH(lptr[2]);
    for (i=0; i<j; i++) second[i] = dptr[i];

    penalty[0] = NUMERIC_VALUE(lptr[3]);

    j = LENGTH(lptr[4]);
    if (j==1) flag[0] = LOGICAL_VALUE(lptr[4]);
    else {
	logical = LOGICAL_POINTER(lptr[4]);
	for (i=0; i<j; i++) flag[i] = logical[i];
	}
    }
