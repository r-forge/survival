/*   SCCS $Id: surv_callback.c,v 1.1 1999-01-31 21:10:49 therneau Exp $
** callback routines for the survreg "other" distributions
** This is modeled directly on the interface code for glm/gam
*/
#include "survS.h"
#include "survproto.h"
#include "eval.h"

static long nframe;
static s_name *survlist;
static s_object *expr1;   /* Code to be executed */

/*
** The first routine just saves away the location of things
*/
void init_survcall(long *ptr1, s_object **ptr2) {
S_EVALUATOR
    long old_ptr;
    s_object *cl; 
   	
    nframe = *ptr1;  /* the frame number of the calling routine */	
	
    /* find the survlist object, and save away it's location */
    survlist = make_name("survlist", S_evaluator);  /*name of desired object*/
    while(nframe > 1 && !(cl= find_in_frame(survlist, nframe, S_evaluator)))
          nframe = parent_frame[nframe];
    if(nframe == 1) PROBLEM
          "survreg C code couldn't find object \"survlist\""
            ERROR;

    /* now grab a copy of the expression */
    old_ptr = set_alloc(nframe, S_evaluator);       /*save prior location */
    expr1 = copy_data(*ptr2, Frames->value.tree[nframe-1], S_evaluator);
    cl = copy_data(cl, NULL_ENTRY, S_evaluator);        /* ? */
    set_in_frame(nframe, cl, survlist, S_evaluator);    /* ? */
    set_alloc(old_ptr, S_evaluator);			/*restore prior */
    }


/*
** This part is called by the survreg5 function, to get the distribution
*/
void surv_callback (double *z, double *dist) {
S_EVALUATOR
    s_object *coxlist, *temp;
    long preva, prevf;

    /* Find the survlist vector, back in the original S calling frame */
    preva = set_alloc(nframe, S_evaluator);

    coxlist = find_in_frame(survlist, nframe, S_evaluator);
    if (!coxlist)  Recover("Couldn't find survlist!", NULL, S_evaluator); 

    /* Now plug the new value of z into it */
    temp = xact_comp(coxlist, "z", S_evaluator);
    if (!temp) Recover("Couldn't find z in the list", NULL, S_evaluator);
    temp = coevec(temp, DOUBLE, TRUE, CHECK_IT, S_evaluator);
    Memmove(temp->value.Double, z, temp->length);

    
   /* evaluate the expression */
    prevf = Nframe;
    set_frame(nframe, S_evaluator);
    eval(expr1, S_evaluator);
    
    /* Grab the updated values from the list */
    
    coxlist = find_in_frame(survlist, nframe, S_evaluator);
    if (!coxlist)  Recover("Couldn't find survlist!", NULL, S_evaluator); 
	
    temp = xact_comp(coxlist, "density", S_evaluator);
    if (!temp) Recover("Couldn't find density in survlist", NULL, S_evaluator);
    temp = coevec(temp, DOUBLE, TRUE, CHECK_IT, S_evaluator);
    Memmove(dist, temp->value.Double, temp->length);

    /* Clean up */
    if (preva > nframe) set_frame(preva, S_evaluator);
    set_alloc(preva, S_evaluator);
    set_frame(prevf, S_evaluator);
    }
