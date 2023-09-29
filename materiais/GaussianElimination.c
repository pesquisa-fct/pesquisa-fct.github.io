#include<stdio.h>
#include<stdlib.h>

void GaussianElimination(double **A, double *b, int n){

	int p, i, j;
	double aux;
	p=0;
	while (p<n-1){
		for (i=p+1; i<n; i++){
			b[i] = b[i] - A[i][p]*(b[p]/A[p][p]);
			aux = A[i][p];
			for (j=p; j<n; j++){
				A[i][j] = A[i][j] - A[p][j]*(aux/A[p][p]);
			}
		}
		p++;
	}
}
