#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "LUdecomposition.h"
#include "SistemasTriangulares.h"
#include "GaussianElimination.h"
#include "GradientMethod.h"
#include "ConjugateGradientMethod.h"
#include "FuncoesAuxiliares.h"
#include "SOR.h"
#include "MetodoPotencia.h"

void main(){

	double **A, *b, *x, detA, *u;
	double **M, **Aq;
	double lambdaMax, lambdaMin, q;
	int n=3, i, j, k, lu, inv;
	int method, methodEigen, iter=0;
	double tol;

	//__________
	//Parameters:
	//__________
	//Tolerancia
	tol = 1.0e-10;
	//__________
	
	//__________
	//Method for linera system:
	//__________
	//LU decomposition	:	0
	//Gaussian Elimination	:	1
	//Cholesky decomposition:	2
	//QR decomposition	:	3
	//SVD decomposition	:	4
	//Jacobi-Richardson	:	5
	//Gauss-Seidel		:	6
	//SOR			:	7
	//Gradientes		:	8
	//Conjugate Gradient	:	9
	method=0;
	//_________

	//_________
	//Method for eigenvalues:
	//_________
	//Metodo da potencia	:	0
	methodEigen = 0;
	//_________	


	//memory allocation
	A =(double **)malloc(sizeof(double *)*n);
	Aq =(double **)malloc(sizeof(double *)*n);
	M = (double **)malloc(sizeof(double *)*(n-1));
	for (k=0; k<n; k++){
		A[k] =(double *)malloc(sizeof(double)*n);
		Aq[k] =(double *)malloc(sizeof(double)*n);
	}	
	for (k=0; k<n-1; k++){
		M[k] = (double *)malloc(sizeof(double)*(n-1));
	}
	b =(double *)malloc(sizeof(double)*n);
	x =(double *)malloc(sizeof(double)*n);
	u =(double *)malloc(sizeof(double)*n);

	
	//For linear systems
	//Matrix definition
	A[0][0] = 10.0;
	A[0][1] = 1.0;
	A[0][2] = 0.0;

	A[1][0] = 1.0;
	A[1][1] = 10.0;
	A[1][2] = 1.0;

	A[2][0] = 0.0;
	A[2][1] = 1.0;
	A[2][2] = 10.0;

	//Independent vector
	b[0] = 11.0;
	b[1] = 11.0;
	b[2] = 1.0;
	

	/*
	//4x4 matrix
	A[0][0] = 2.0;
	A[0][1] = 3.0;
	A[0][2] = - 1.0;
	A[0][3] = 2.0;
	      				
	A[1][0] = 0.0;
	A[1][1] = 4.0;
	A[1][2] = -3.0;
	A[1][3] = 5.0;
	
	A[2][0] = 1.0;
	A[2][1] = 2.0;
	A[2][2] = 1.0;
	A[2][3] = 3.0;

	A[3][0] = 0.0;
	A[3][1] = 4.0;
	A[3][2] = 1.0;
	A[3][3] = 0.0;
	*/
	/*
	//5x5 matrix
	A[0][0] = 2.0;
	A[0][1] = 3.0;
	A[0][2] = - 1.0;
	A[0][3] = 2.0;
	A[0][4] = 2.0;
	      				
	A[1][0] = 0.0;
	A[1][1] = 4.0;
	A[1][2] = -3.0;
	A[1][3] = 5.0;
	A[1][4] = 5.0;
	
	A[2][0] = 1.0;
	A[2][1] = 2.0;
	A[2][2] = 1.0;
	A[2][3] = 3.0;
	A[2][4] = 5.0;

	A[3][0] = 0.0;
	A[3][1] = 4.0;
	A[3][2] = 1.0;
	A[3][3] = 0.0;
	A[3][4] = 5.0;

	A[4][0] = 1.0;
	A[4][1] = 2.0;
	A[4][2] = 3.0;
	A[4][3] = 4.0;
	A[4][4] = 5.0;
	*/

	
	//matrix for eigenvalues testes:
	A[0][0] = 3.0;
	A[0][1] = 0.0;
	A[0][2] = 1.0;
	
	A[1][0] = 2.0;
	A[1][1] = 2.0;
	A[1][2] = 2.0;

	A[2][0] = 4.0;
	A[2][1] = 2.0;
	A[2][2] = 5.0;
	

	/*
	A[0][0] = 3.0;
	A[0][1] = 1.0;
	A[0][2] = 0.0;
	
	A[1][0] = 1.0;
	A[1][1] = 2.0;
	A[1][2] = -1.0;

	A[2][0] = 0.0;
	A[2][1] = -1.0;
	A[2][2] = 0.0;
	*/

	/*
	A[0][0] = 2.0;
	A[0][1] = 1.0;
	A[0][2] = 0.0;
	
	A[1][0] = 2.0;
	A[1][1] = 5.0;
	A[1][2] = 3.0;

	A[2][0] = 0.0;
	A[2][1] = 1.0;
	A[2][2] = 6.0;
	*/

	/*
	n=2;
	A[0][0] = 3.0;
	A[0][1] = 4.0;
	A[1][0] = 2.0;
	A[1][1] = 1.0;
	*/

	//______________________________________
	//Potencia - Maior autovalor		:	0
	//Potencia Inversa - menor autovalor	:	1
	inv = 0;
	//______________________________________
	lambdaMax = MetodoPotencia(A, u, n, tol, inv);
	printf("\nAutovetor\n");
	for (i=0 ;i<n; i++){
		printf("u[%d]=%f\n", i, u[i]);
	}
	printf("\nAutovalor Max\n");
	printf("lambdaMax=%f\n", lambdaMax);//getchar();
	//Deslocamento:
	q = lambdaMax;
	MatrizDeslocada(A, Aq, n, q);
	inv = 0;
	lambdaMin = MetodoPotencia(Aq, u, n, tol, inv) + q;
	printf("\nAutovalor Min\n");
	printf("lambdaMin=%f\n", lambdaMin);//getchar();

	//exit(1);
	
	detA = Determinant(A, n);
	printf("detA=%f\n", detA);//getchar();

	//chute inicial
	for (i=0; i<n; i++){
		x[i] = 0.0;
	}

	switch(method){
		case 0://LU decomposition
			printf("\nLU decomposition\n");
			lu=1;//Esse parametro considera o sistema triangular Inferior com 1's na diagonal principal.
			detA = LUdecomposition(A, n);
			//Solucao dos sistemas
			SistemaTriangularInferior(A, b, x, n, lu);
			SistemaTriangularSuperior(A, x, x, n);//O vetor independente eh a solucao do sistema triangular inferior.
			break;
		case 1://Gaussian Elimination
			printf("\nGaussian Elimination\n");
			GaussianElimination(A, b, n);
			SistemaTriangularSuperior(A, b, x, n);
			break;
		case 2://Cholesk Decomposition
			break;
		case 3://QR decomposition
			break;
		case 4://SVD decomposition
			break;
		case 5://Jacobi-Richardson
			break;
		case 6://Gauss-Seidel
			break;
		case 7://SOR
			printf("\nSOR method\n");
			iter = SOR(A, b, x, n, tol);
			break;
		case 8://Gradientes
			printf("\nGradient Method\n");
			iter = GradientMethod(A, b, x, n, tol);
			break;
		case 9://Gradientes Conjugados
			printf("\nConjugate gradient method\n");
			iter = ConjugateGradientMethod(A, b, x, n, tol);
			break;
		default:
			break;

	}

	//Imprimir solucao final
	printf("\nSolucao final:\n");
	for (k=0; k<n; k++){
		printf("x[%d]=%f\n", k, x[k]);
	}
	//Numero de iteracoes
	printf("\ndetA=%lf\n", detA);
	printf("\nNumero de iteracoes: iter=%d\n", iter);

	//free memory
	for (k=0; k<n; k++){
		free(A[k]);
		free(Aq[k]);
	}
	for (k=0; k<n-1; k++){
		free(M[k]);
	}
	free(A);
	free(Aq);
	free(M);
	free(b);
	free(x);
	free(u);

}
