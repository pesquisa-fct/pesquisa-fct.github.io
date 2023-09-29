#include "GradientMethod.h"
#include "FuncoesAuxiliares.h"

int GradientMethod(double **A, double *b, double *v, int n, double tol){

	int i, j, iter=0;;
	double r[n], vOld[n], tmin, residuo, erro, prod[n], dif[n], sum;

	//Inicial condition (chute inicial)
	//for (i=0; i<n; i++){
		//v[i] = 0.0;
	//}

	//Calculo do residuo inicial:
	for (i=0; i<n; i++){
		sum = 0.0;
		for (j=0; j<n; j++){
			sum += A[i][j]*v[j];
		}
		r[i] = sum - b[i];
	}

	residuo = Norm(r, n);
	erro = 10.0;



	while (residuo > tol & erro > tol){

		iter++;
		//Copia solucao nova para velha
		for (i=0; i<n; i++){
			vOld[i] = v[i];
		}
		//Matrix product A*r
		for (i=0; i<n; i++){
			sum = 0.0;
			for (j=0; j<n; j++){
				sum+=A[i][j]*r[j];
			}
			prod[i] = sum;
		}

		//calculo do t_min
		tmin = InnerProduct(r,r,n)/InnerProduct(prod, r, n);

		//Calculo da nova aproximacao e do novo residuo
		for (i=0; i<n; i++){
			v[i] = v[i] - tmin*r[i];
			r[i] = r[i] - tmin*prod[i];
			dif[i] = v[i] - vOld[i];
		}

		residuo = Norm(r, n);
		erro = Norm(dif, n)/Norm(v, n);
	}

	return iter;
}
