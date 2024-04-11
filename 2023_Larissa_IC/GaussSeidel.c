#include "GaussSeidel.h"

int GaussSeidel(double **A, double *b, double *x, int n){
	int i, j, iter;
	double aux, tol = 1.0e-04, RelErr, xOld[n], soma, erro[n];

	//Dividir pela diagonal principal (A^* e b^*)
	for (i=0; i<n; i++){
		aux = A[i][i];
		for (j=0; j<n; j++){
			A[i][j] = A[i][j]/aux;
		}
		b[i] = b[i]/aux;
	}

	//Algoritmo
	iter = 1;
	RelErr = 10.0;
	while (RelErr > tol){

		//copia solucao nova para velha
		for (i=0; i<n; i++){
			xOld[i] = x[i];
		}

		//algoritmo
		for (i=0; i<n; i++){
			soma = 0.0;
			for (j=0; j<n; j++){
				if (i != j){
					soma += A[i][j]*x[j];
				}
			}
			x[i] = b[i] - soma;
			erro[i] = x[i] - xOld[i];
		}

		RelErr = Norm(erro, n)/Norm(x, n);
		iter++;
	}


	return (iter);
}
