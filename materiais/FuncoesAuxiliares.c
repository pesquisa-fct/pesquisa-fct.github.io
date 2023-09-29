#include "FuncoesAuxiliares.h"

double Norm(double *v, int n){

	double value=0.0;

	for (int i=0; i<n; i++){
		value += value + v[i]*v[i];
	}
	value = sqrt(value);

	return (value);
}

double MaxNorm(double *v, int n){

	int i;
	double max, value;
	max = fabs(v[0]);
	for (i=1; i<n; i++){
		value = fabs(v[i]);
		if (value > max){
			max = value;
		}
	}

	return (max);
}

void MatrizDeslocada(double **A, double **Aq, int n, double q){
	int i, j;
	for (i=0; i<n; i++){
		for (j=0; j<n; j++){
			Aq[i][j] = A[i][j];
			if (i==j){
				Aq[i][j] += -q;
			}
		}
	}
}

double InnerProduct(double *u, double *v, int n){

	double value=0.0;
	int i;

	for (i=0; i<n; i++){
		value += u[i]*v[i];
	}

	return value;
}

double *MatrixProduct(double **A, double *v, int n){
	int i, j;
	double sum, *prod;
	prod = (double *)malloc(sizeof(double)*n);

	for (i=0; i<n; i++){
		sum = 0.0;
		for (j=0; j<n; j++){
			sum += A[i][j]*v[j];
		}
		prod[i] = sum;
	}

	return prod;

	free(prod);
}

//Essa funcao tem por objetivo calcular as matrizes menores, usadas no calculo do determinante, por meio 
//do teorema de Laplace. Tambem, usada no calcula da matriz dos cofatores, adjunta e inversa de A.
void MenorMatriz(double **A, double **M, int i, int j, int n){
	int ii, jj, dim;

	//Parametros
	//___________
	//Dimensao da matriz menor M
	dim = n-1;

	//Bloco I
	for (ii=0; ii<i; ii++){
		for (jj=0; jj<j; jj++){
			M[ii][jj] = A[ii][jj];
		}
	}

	//Bloco II
	for (ii=0; ii<i; ii++){
		for (jj=j; jj<dim; jj++){
			M[ii][jj] = A[ii][jj+1];
		}
	}
	
	//Bloco III
	for (ii=i; ii<dim; ii++){
		for (jj=0; jj<j; jj++){
			M[ii][jj] = A[ii+1][jj];
		}
	}

	//Bloco IV
	for (ii=i; ii<dim; ii++){
		for (jj=j; jj<dim; jj++){
			M[ii][jj] = A[ii+1][jj+1];
		}
	}


}

double Determinant(double **A, int n){

	double detA;
	int i, j;

	if (n == 1){
		detA = A[0][0];
	}else if (n==2){
		detA = A[0][0]*A[1][1] - A[0][1]*A[1][0];
	}else if (n==3){
		detA = A[0][0]*A[1][1]*A[2][2] + A[0][1]*A[1][2]*A[2][0] + A[0][2]*A[1][0]*A[2][1] - (A[0][2]*A[1][1]*A[2][0] + A[0][1]*A[1][0]*A[2][2] + A[0][0]*A[1][2]*A[2][1]);
	}else if (n>3){
		//esolher linha
		//____________
		i = n-1;
		//____________

		double **M;
		int dim=n-1;
		//memory allocation
		M = (double **)malloc(sizeof(double *)*dim);
		for (i=0; i<dim; i++){
			M[i] = (double *)malloc(sizeof(double)*dim);
		}

		detA = 0.0;
		for (j=0; j<n; j++){
			MenorMatriz(A, M, i, j, n);
			if ( ( (i+j)%2 ) == 0){
				detA += A[i][j]*Determinant(M, dim);
			}else{
				detA += - A[i][j]*Determinant(M, dim);
			}
		}

		//free memory
		for (i=0; i<dim; i++){
			free(M[i]);
		}
		free(M);
	}


	return (detA);

}
