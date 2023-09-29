#include "SOR.h"

int SOR(double **A, double *b, double *x, int n, double tol){

	int iter, i, j;
	double soma, w, erro, vec[n], xOld[n], aux;

	w=1.2;//Parametro de relaxacao

	//inicial solution
	//for (i=0; i<n ;i++){
		//x[i] = 0.0;
	//}

	//Dividir a matriz dos coeficientes e b pela diagonal principal
	for (i=0; i<n; i++){
		aux = A[i][i];
		for (j=0; j<n; j++){
			A[i][j] = A[i][j]/aux;
		}
		b[i] = b[i]/aux;
	}

	erro=10;
	//Algoritmo SOR
	while (erro > tol){
		//copia solucao nova para velha
		for (i=0; i<n; i++){
			xOld[i] = x[i];
		}
		//Metodo SOR
		for (i=0; i<n; i++){
			soma =0.0;
			for (j=0; j<n; j++){
				soma += A[i][j]*x[j];
			}
			x[i] = x[i] + w*(-soma + b[i]);
			vec[i] = x[i] - xOld[i];

		}

		erro = Norm(vec, n)/Norm(x, n);
		iter++;
	}		


	return iter;
}
