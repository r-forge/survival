/* SCCS: $Id: chsolve2.c,v 1.2 1991-06-12 14:00:19 therneau Exp $ */
/*
** Solve the equation Ab = y, where the cholesky decomposition of A and y
**   are the inputs.
**
** Input  **matrix, which contains the chol decomp of an n by n
**   matrix in its lower triangle.
**        y[n] contains the right hand side
**
**  y is overwriten with b
**
**  Terry Therneau
*/

chsolve (matrix, n, y)
int  n;
double **matrix, y[];
     {
     register int i,j;
     register double temp;

     /*
     ** solve Fb =y
     */
     for (i=0; i<n; i++) {
	  temp = y[i] ;
	  for (j=0; j<i; j++)
	       temp -= y[j] * matrix[i][j] ;
	  y[i] = temp ;
	  }
     /*
     ** solve DF'z =b
     */
     for (i=(n-1); i>=0; i--) {
	  temp = y[i]/matrix[i][i];
	  for (j= i+1; j<n; j++)
	       temp -= y[j]*matrix[j][i];
	  y[i] = temp;
	  }
     }
