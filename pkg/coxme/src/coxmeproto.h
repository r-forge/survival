/*
** This contains the prototype calls for all the .c functions that
**    are called by another C function
**  Mostly, it stops errors due to having things declared differently
**    in different routines.
*/

/* 
 * This redefine is a temporary fix.  Since Splus7 added the kinship library
 * into the Splus core, calling cholesky4 from inside other C routines
 * called Splus's version and not the kinship library's.  This redefine 
 * bypasses Splus's cholesky4 by renaming kinship's version.
 * Eric Lunde, 2006-04-04
*/
#define cholesky4 cholesky4_fix_for_splus7

void bdsmatrix_prod2(int nblock,     int *bsize,     int nrow,
		     double *bmat,   double *rmat,  
		     double *y,      double *result, int *itemp) ;

void chinv4(double **matrix, 	int n, 	 	int nblock, 	int *bsize, 
	    double *bd,	        int flag) ;


void chinv5(double **matrix ,	int n,		int flag);

int cholesky4(double **matrix,	int n,		int nblock, 	int *bsize,
	      double *bd, 	double toler) ;

int cholesky5(double **matrix, 	int n, 		double toler);

void chsolve4(double **matrix, 	int n, 		int nblock, 	int *bsize,
	      double *bd, 	double *y, 	int flag);

void chsolve5(double **matrix, 	int n, 		double *y, 	int flag);

double **dmatrix(double *array, int ncol, 	int nrow);

void bdsmatrix_prod4(int nrow,    int nblock,   int *bsize, 
		    double *bmat, double *rmat,    
		    int nfrail,   double *y);




