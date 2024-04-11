#include "ConjugateGradientMethod.h"

int ConjugateGradientMethod(double **A, double *b, double *v, int n, double tol){

	int i, j, iter;
	double r[n], p[n], dif[n], vOld[n], alpha, q, *prod;
	double sum, residuo, erro;

	prod = (double *)malloc(sizeof(double)*n);

	//Inicial cicle
	//____________________
	//Inicial approximation v_0
	//for (i=0; i<n; i++){
		//v[i] = 0.0;
	//}
	//____________________
	//Residuo
	for (i=0; i<n; i++){
		sum = 0.0;
		for (j=0; j<n; j++){
			sum += A[i][j]*v[j];
		}
		r[i] = sum - b[i];
		//printf("A*v[%d]=%f\n", i, sum);getchar();
	}
	//____________________
	//Inicial Direction
	for (i=0; i<n; i++){
		p[i] = -r[i];
	}
	//Produto matricial
	for (i=0; i<n; i++){
		sum = 0.0;
		for (j=0; j<n; j++){
			sum += A[i][j]*p[j];
		}
		prod[i] = sum;
	}
	alpha = - InnerProduct(r, p, n)/InnerProduct(prod, p, n);
	for (i=0; i<n; i++){
		v[i] = v[i] + alpha*p[i];
	}
	//Atualizando residuo
	for (i=0; i<n; i++){
		sum =0.0;
		for (j=0; j<n; j++){
			sum += A[i][j]*v[j];
		}
		r[i] = sum - b[i];
	}

	iter = 1;
	residuo = Norm(r, n);
	erro = 10.0;

	while (residuo > tol && erro > tol){

		//Copia solucao nova para velha:
		for (i=0; i<n; i++){
			vOld[i] = v[i];
		}

		//Auxiliar alpha
		alpha = InnerProduct(MatrixProduct(A, r, n), p, n)/InnerProduct(MatrixProduct(A, p, n), p, n);

		//New direction p
		for (i=0; i<n; i++){
			p[i] = -r[i] + alpha*p[i];
		}

		//Auxiliar q
		q = -InnerProduct(r, p, n)/InnerProduct(MatrixProduct(A, p, n), p, n);

		//New approximation
		for (i=0; i<n; i++){
			v[i] = v[i] + q*p[i];
		}

		//New residuo
		prod = MatrixProduct(A, p, n);
		for (i=0; i<n; i++){
			r[i] = r[i] + q*prod[i];
			dif[i] = v[i] - vOld[i];
		}

		residuo = Norm(r, n);
		erro = Norm(dif, n)/Norm(v, n);
		iter++;
		
	}

	free(prod);

	return iter;
}
