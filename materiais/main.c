#include<stdio.h>
#include<stdlib.h>
#include<math.h>
//#include "InnerProduct.h"
#include "GaussianElimination.h"
#include "BackwardSubstitution.h"
//#include "LU_Decomposition.h"
#include "JacobiRichardson.h"

int main(){
	int i, j;
	int n=3, method;
	double **A, **L, **U; 
	double *x, *b;

	//Memory allocation
	A = (double **)malloc(sizeof(double *)*(n+1));
	L = (double **)malloc(sizeof(double *)*(n+1));
	U = (double **)malloc(sizeof(double *)*(n+1));
	x = (double *)malloc(sizeof(double)*(n+1));
	b = (double *)malloc(sizeof(double)*(n+1));

	for (i=0; i<n+1; i++){
		A[i] = (double *)malloc(sizeof(double)*(n+1));
		L[i] = (double *)malloc(sizeof(double)*(n+1));
		U[i] = (double *)malloc(sizeof(double)*(n+1));
	}

	//Definir matriz A
	A[0][0] = 10.0;
	A[0][1] = 2.0;
	A[0][2] = 1.0;
	//A[0][3] = 1.0;
	b[0] = 7.0;

	A[1][0] = 1.0;
	A[1][1] = 5.0;
	A[1][2] = 1.0;
	//A[1][3] = 0.0;
	b[1] = -8.0;

	A[2][0] = 2.0;
	A[2][1] = 3.0;
	A[2][2] = 10.0;
	//A[2][3] = 0.0;
	b[2] = 6.0;

	//A[3][0] = 0.0;
	//A[3][1] = 0.0;
	//A[3][2] = 1.0;
	//A[3][3] = 1.0;
	//b[3] = 700.0;

	//_______________________
	//method:
	//Gaussian Elimination:	0
	//Jacobi-Richardson:	1
	method = 1;
	//-----------------------

	switch(method){
		case 0: //Gaussian Elimination
			printf("\nMethod: Gaussian Elimination\n");
			GaussianElimination(A, U, b, n);//U é a matriz transformada
			BackwardSubstitution(A, x, b, n);	
		break;

		case 1: //Jacobi-Richardson
			printf("\nMethod: Jacobi-Richardson\n");
			JacobiRichardson(A, x, b, n);
			break;
	}

	//Print the solution
	for (i=0; i<n; i++){
		printf("x[%d]=%f \n", i, x[i]);
	}

	//Free memory
	for (i=0; i<(n+1); i++){
		free(A[i]);
		free(L[i]);
		free(U[i]);
	}
	free(A);
	free(L);
	free(U);
	free(x);

	return (0);
}
