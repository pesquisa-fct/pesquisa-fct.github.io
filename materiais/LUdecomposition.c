#include<stdio.h>
#include<stdlib.h>

double LUdecomposition(double **A, int n){

	int i, j, k;
	double sum, detA1, detA2, detA;	

	//Verificando as hipoteses do teorema da decomposicao LU para n<=3.
	if (n <= 3){
		if (n==3){
			detA1 = A[0][0];
			detA2 = A[0][0]*A[1][1] - A[0][1]*A[1][0];
		}else if (n==2){
			detA1 = A[0][0];
			detA2 = 10.0;
		}else{
			printf("n=1; caso trivial\n");
			exit(1);
		}

		if ((detA1 == 0.0) && (detA2 == 0.0)){
			printf("Um dos menores principais nao eh diferente de zero\n");
			exit(1);
		}
		else{
			printf("Eh possivel aplicar a decomposicao LU\n");
		}
	}

	//Decomposicao LU
	for (i=0; i<n; i++){
		for (j=0; j<n; j++){

			//Matrix U
			if ( (i < j) || (i==j)){
				sum = 0.0;
				if (i != 0){
					for (k=0; k<i; k++){
						sum = sum + A[i][k]*A[k][j];
					}
				}
				A[i][j] = A[i][j] - sum;//U[i][j] = A[i][j] - sum_{k=1}^{i-1} L[i][k]*U[k][j];
			}
			//Matrix L
			else if (i > j){
				sum = 0.0;
				if (j>0){
					for (k=0; k<j; k++){
						sum = sum + A[i][k]*A[k][j];

					}
				}
				A[i][j] = (A[i][j] - sum)/A[j][j];
			}

		}//end j
	}//end i

	//Calculo do determinante a partir da matriz U
	detA = 1.0;
	for (k=0; k<n; k++){
		detA = detA*A[k][k];
	}

	return detA;

}
