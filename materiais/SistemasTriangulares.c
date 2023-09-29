#include<stdio.h>
#include<stdlib.h>

void SistemaTriangularInferior(double **A, double *b, double *x, int n, int lu){

	int i, k;
	double sum;

	//Matriz triangular inferior com diagonal igual a 1
	if (lu != 0){
		x[0] = b[0];
		for (i=1; i<n; i++){
			sum = 0.0;
			for (k=0; k<i; k++){
				sum += A[i][k]*x[k];
			}
			x[i] = b[i] - sum;
		}
	}//Matriz triangular inferior convencional
	else{
		x[0] = b[0]/A[0][0];
		for (i=1; i<n; i++){
			sum = 0.0;
			for (k=0; k<i; k++){
				sum += A[i][k]*x[k];
			}
			x[i] = (b[i] - sum)/A[i][i];
		}

	}

	
}

void SistemaTriangularSuperior(double **A, double *b, double *x, int n){

	int i, k;
	double sum;

	x[n-1] = b[n-1]/A[n-1][n-1];
	for (i=2; i<=n; i++){
		sum=0.0;
		for (k=n-i+1; k<n; k++){
			sum += A[n-i][k]*x[k];
		}
		x[n-i] = (b[n-i] - sum)/A[n-i][n-i];
	}

}
