/*
** callback routines for the coxph frailty interface
** This is modeled directly on the interface code for glm/gam
*/
#include "survS.h"
#include "survproto.h"
#include "eval.h"

static long nframe;
static s_name *coxlist1, *coxlist2;
static s_object *expr1;   /* penalized term(s) */
static s_object *expr2;   /* sparse term */

/*
** The first routine just saves away the location of things
*/
void init_coxcall1(long *ptr1, s_object **ptr2) {
S_EVALUATOR
    long old_ptr;
    s_object *cl; 
   	
    nframe = *ptr1;  /* the frame number of the calling routine */	
	
    /* find the coxlist1 object, and save away it's location */
    coxlist1 = make_name("coxlist1", S_evaluator);  /*name of desired object*/
    while(nframe > 1 && !(cl= find_in_frame(coxlist1, nframe, S_evaluator)))
          nframe = parent_frame[nframe];
    if(nframe == 1) PROBLEM
          "coxph C code couldn't find object \"coxlist1\""
            ERROR;

    /* now grab a copy of the expression */
    old_ptr = set_alloc(nframe, S_evaluator);       /*save prior location */
    expr1 = copy_data(*ptr2, Frames->value.tree[nframe-1], S_evaluator);
    cl = copy_data(cl, NULL_ENTRY, S_evaluator);        /* ? */
    set_in_frame(nframe, cl, coxlist1, S_evaluator);    /* ? */
    set_alloc(old_ptr, S_evaluator);			/*restore prior */
}

void init_coxcall2(long *ptr1, s_object **ptr2) {
S_EVALUATOR
    long old_ptr;
    s_object *cl; 
   	
    nframe = *ptr1;  /* the frame number of the calling routine */	
	
    coxlist2 = make_name("coxlist2", S_evaluator);
    while(nframe > 1 && !(cl =find_in_frame(coxlist2, nframe, S_evaluator)))
          nframe = parent_frame[nframe];
    if(nframe == 1) PROBLEM
          "coxph C code couldn't find object \"coxlist2\""
            ERROR;

     /* now grab a copy of the expression */
    old_ptr = set_alloc(nframe, S_evaluator);  /*save prior location */
    expr2 = copy_data(*ptr2, Frames->value.tree[nframe-1], S_evaluator);/*get*/
    cl = copy_data(cl, NULL_ENTRY, S_evaluator);        /* ? */
    set_in_frame(nframe, cl, coxlist2, S_evaluator);    /* ? */
    set_alloc(old_ptr, S_evaluator);  /*restore prior */
}

/*
** This part is called by the coxfit4 function, to get the penalty terms
*/
void cox_callback (int which, double *coef, double *first, 
	           double *second, double *penalty, long *flag) {
S_EVALUATOR
    s_object *coxlist, *temp;
    long preva, prevf;

    /* Find the coxlist vector, back in the original S calling frame */
    preva = set_alloc(nframe, S_evaluator);

    if (which==1) {
        coxlist = find_in_frame(coxlist1, nframe, S_evaluator);
	if (!coxlist)  Recover("Couldn't find coxlist1!", NULL, S_evaluator); 
	}
    else {
        coxlist = find_in_frame(coxlist2, nframe, S_evaluator);
	if (!coxlist)  Recover("Couldn't find coxlist2!", NULL, S_evaluator); 
        }

    /* Now plug the new value of coef into it */
    temp = xact_comp(coxlist, "coef", S_evaluator);
    if (!temp) Recover("Couldn't find coef in the list", NULL, S_evaluator);
    temp = coevec(temp, DOUBLE, TRUE, CHECK_IT, S_evaluator);
    Memmove(temp->value.Double, coef, temp->length);

    
   /* evaluate the expression */
    prevf = Nframe;
    set_frame(nframe, S_evaluator);
    if (which==1) eval(expr1, S_evaluator);
    else          eval(expr2, S_evaluator);

    /* Grab the updated values from the list */
    if (which==1) {
        coxlist = find_in_frame(coxlist1, nframe, S_evaluator);
	if (!coxlist)  Recover("Couldn't find coxlist1!", NULL, S_evaluator); 
	}
    else {
        coxlist = find_in_frame(coxlist2, nframe, S_evaluator);
	if (!coxlist)  Recover("Couldn't find coxlist2!", NULL, S_evaluator); 
        }
	
    temp = xact_comp(coxlist, "coef", S_evaluator);
    if (!temp) Recover("Couldn't find coef in coxlist", NULL, S_evaluator);
    temp = coevec(temp, DOUBLE, TRUE, CHECK_IT, S_evaluator);
    Memmove(coef, temp->value.Double, temp->length);
    
    temp = xact_comp(coxlist, "first", S_evaluator);
    if (!temp) Recover("Couldn't find first in coxlist", NULL, S_evaluator);
    temp = coevec(temp, DOUBLE, TRUE, CHECK_IT, S_evaluator);
    Memmove(first, temp->value.Double, temp->length);

    temp = xact_comp(coxlist, "second", S_evaluator);
    if (!temp) Recover("Couldn't find second in coxlist", NULL, S_evaluator);
    temp = coevec(temp, DOUBLE, TRUE, CHECK_IT, S_evaluator);
    Memmove(second, temp->value.Double, temp->length);
    
    temp = xact_comp(coxlist, "flag", S_evaluator);
    if (!temp) Recover("Couldn't find flag in coxlist", NULL, S_evaluator);
    temp = coevec(temp, LGL, TRUE, CHECK_IT, S_evaluator);
    Memmove(flag, temp->value.Long, temp->length);

    temp = xact_comp(coxlist, "penalty", S_evaluator);
    if (!temp) Recover("Couldn't find penalty in coxlist", NULL, S_evaluator);
    temp = coevec(temp, DOUBLE, TRUE, CHECK_IT, S_evaluator);
    Memmove(penalty, temp->value.Double, temp->length);

    /* Clean up */
    if (preva > nframe) set_frame(preva, S_evaluator);
    set_alloc(preva, S_evaluator);
    set_frame(prevf, S_evaluator);
    }
