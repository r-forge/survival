/* SCCS $Id: cholesky2.c,v 4.1 1992-03-04 16:51:48 therneau Exp $  */
/*
** subroutine to do Cholesky decompostion on a matrix: C = FDF'
**   where F is lower triangular with 1's on the diagonal, and D is diagonal
**
** arguments are:
**     n         the size of the matrix to be factored
**     **matrix  a ragged array containing an n by n submatrix to be factored
**
**  The factorization is returned in the lower triangle, D occupies the
**    diagonal and the upper triangle is left undisturbed.  The original
**    contents of the lower triangle are ignored, and could be zeros.
**
**  0 is returned upon successful factorization, the index of the
**    offending column is returned upon failure
**
**   Terry Therneau
*/

cholesky(matrix, n)
double **matrix;
int n;
     {
     register double temp;
     register int  i,j,k;

     for (i=0; i<n; i++) {
	  for (j=0; j<i; j++) {
	       temp =0;
	       for (k=0; k<j; k++)
		    temp += matrix[j][k] *matrix[i][k]*matrix[k][k]  ;
	       matrix[i][j] = (matrix[j][i] -temp) / matrix[j][j];
	       }

	  temp =0;
	  for (k=0; k<i; k++)
	       temp += matrix[i][k]*matrix[i][k]*matrix[k][k];
	  matrix[i][i] -= temp;
	  if ( matrix[i][i] <= 0.0)   return(i+1);
	  }
     return(0);
     }
