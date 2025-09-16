#include<stdio.h>
#include<stdlib.h>
#include<math.h>

void DecomposicaoQR(double **A, double **Q, double **R, int n);
void MetodoQR(double **A, double **R, double **Q, double tol, int n);

int main(){

	double **A, **Q, **R, soma, tol=1.0e-04;
	int n=3, i, j, k;

	A = (double **)malloc(sizeof(double *)*n);
	Q = (double **)malloc(sizeof(double *)*n);
	R = (double **)malloc(sizeof(double *)*n);

	for (i=0; i<n; i++){
		A[i] = (double *)malloc(sizeof(double)*n);
		Q[i] = (double *)malloc(sizeof(double)*n);
		R[i] = (double *)malloc(sizeof(double)*n);
	}

	A[0][0] = 2.0;
	A[0][1] = 1.0;
	A[0][2] = 0.0;

	A[1][0] = 2.0;
	A[1][1] = 5.0;
	A[1][2] = 3.0;

	A[2][0] = 0.0;
	A[2][1] = 1.0;
	A[2][2] = 6.0;

	//DecomposicaoQR(A, Q, R, n);
	MetodoQR(A, R, Q, tol, n);
	/*
	for (i=0; i<n; i++){
		for (j=0; j<n; j++){
			printf("Q[%d][%d]=%f\t", i,j,Q[i][j]);
		}
		printf("\n");
	}
	for (i=0; i<n; i++){
		for (j=0; j<n; j++){
			printf("R[%d][%d]=%f\t", i,j,R[i][j]);
		}
		printf("\n");
	}
	*/

	/*
	//Reconstruindo a matriz A:
	for (i=0; i<n; i++){
		for (j=0; j<n; j++){
			soma = 0.0;
			for (k=0; k<n; k++){
				soma = soma + Q[i][k]*R[k][j];
			}
			A[i][j] = soma;
			printf("A[%d][%d]=%f\t", i, j, A[i][j]);
		}
		printf("\n");
	}
	*/


	for (i=0; i<n; i++){
		free(A[i]);
		free(Q[i]);
		free(R[i]);
	}
	free(A);
	free(Q);
	free(R);



	return (0);
}

void DecomposicaoQR(double **A, double **Q, double **R, int n){
	int i, j, p, q, k;
	double **U, **P;
	double c, s, r, soma, aux;

	U = (double **)malloc(sizeof(double *)*n);
	P = (double **)malloc(sizeof(double *)*n);


	for (i=0; i<n; i++){
		U[i] = (double *)malloc(sizeof(double)*n);
		P[i] = (double *)malloc(sizeof(double)*n);

	}

	//Inicializar U:
	for (i=0; i<n; i++){
		for (j=0; j<n; j++){
			R[i][j] = A[i][j];
			//printf("Rij = %f\n", R[i][j]);getchar();
			if (i==j){
				P[i][j] = Q[i][j] = U[i][j] = 1.0;
			}else{
				P[i][j] = Q[i][j] = U[i][j] = 0.0;
			}
		}
	}

	for (p=0; p<(n-1); p++){
		for (q=p+1; q<n; q++){
			if (R[q][p] != 0.0){
				//printf("matriz de rotacao U != I\n");
				//printf("p=%d, q=%d\n", p, q);

				r = sqrt( R[q][p]*R[q][p] + R[p][p]*R[p][p] );
				c = R[p][p]/r;
				s = R[q][p]/r;
				U[p][p] = c;
				U[p][q] = s;
				U[q][p] = -s;
				U[q][q] = c;
			}else{
				//printf("matriz de rotacao U=I\n");
				//printf("p=%d, q=%d\n", p, q);
			}
			/*
			for (i=0;i<n;i++){
				for (j=0;j<n;j++){
					printf("U[%d][%d]=%f\t", i,j,U[i][j]);
				}
				printf("\n");
			}getchar();
			*/

			/*
			printf("Antes Matriz R:\n");
			for (i=0; i<n; i++){
				for (j=0; j<n; j++){
					printf("R[%d][%d]=%f\t", i, j, R[i][j]);
				}
				printf("\n");
			}//getchar();
			*/

			for (j=0; j<n; j++){

				//Estou usando o P como variavel auxiliar, para nao estragar a matriz R.
				//Linha p:
				k=p;
				//aux = R[k][j];
				P[p][j] = U[p][k] * R[k][j];
				k=q;
				P[p][j] += U[p][k] * R[k][j];

				//Linha q
				k=p;
				P[q][j] = U[q][k] * R[k][j];
				k=q;
				P[q][j] += U[q][k] * R[k][j];

				//Atribuo os valores calculados de P para R:
				//Linha p:
				R[p][j] = P[p][j];
				//Linha q
				R[q][j] = P[q][j];

			}

			/*
			printf("Depois Matriz R:\n");
			for (i=0; i<n; i++){
				for (j=0; j<n; j++){
					printf("R[%d][%d]=%f\t", i, j, R[i][j]);
				}
				printf("\n");
			}getchar();
			*/


			//Calculo da Q:
			//printf("Matriz Q\n");
			for (i=0; i<n; i++){
				for (j=0; j<n; j++){
					soma=0.0;
					for (k=0; k<n; k++){
						//printf("Qik=%f, Ujk=%f\n", Q[i][k], U[j][k]);getchar();
						soma = soma + Q[i][k]*U[j][k];
						//printf("soma=%f\n", soma);
					}
					P[i][j] = soma;
					//printf("P[%d][%d]=%f\t", i, j, P[i][j]);

				}
				//printf("\n");
			}//getchar();

			for (i=0; i<n; i++){
				for (j=0; j<n; j++){
					Q[i][j] = P[i][j];
				}
			}

			//Restart na matriz U
			U[p][p] = 1.0;
			U[p][q] = 0.0;
			U[q][p] = 0.0;
			U[q][q] = 1.0;
		}

	}

	for (i=0; i<n; i++){
		free(U[i]);
		free(P[i]);
	}
	free(U);
	free(P);

}

void MetodoQR(double **A, double **R, double **Q, double tol, int n){

	double max=2.0*tol, **U, **P, soma;
	int i, j, k;

	U = (double **)malloc(sizeof(double *)*n);
	P = (double **)malloc(sizeof(double *)*n);

	for (i=0; i<n; i++){
		U[i] = (double *)malloc(sizeof(double)*n);
		P[i] = (double *)malloc(sizeof(double)*n);

	}

	while (max > tol){

		DecomposicaoQR(A, Q, R, n);

		//R*Q
		for (i=0; i<n; i++){
			for (j=0; j<n; j++){
				soma=0.0;
				for (k=0; k<n; k++){
					soma = soma + R[i][k]*Q[k][j];
				}
				A[i][j] = soma;
			}
		}

		/*
		//A*Q
		for (i=0; i<n; i++){
			for (j=0; j<n; j++){
				soma=0.0;
				for (k=0; k<n; k++){
					soma = soma + A[i][k]*Q[k][j];
				}
				U[i][j] = soma;
			}
		}

		//Q^T*(A*Q):
		for (i=0; i<n; i++){
			for (j=0; j<n; j++){
				soma=0.0;
				for (k=0; k<n; k++){
					soma = soma + Q[k][i]*U[k][j];
				}
				A[i][j] = soma;
			}
		}
		*/

		//imprimindo A:
		for (i=0; i<n; i++){
			for (j=0; j<n; j++){
				printf("A[%d][%d]=%f\t", i,j,A[i][j]);
			}
			printf("\n");
		}getchar();

		//Determinando max:
		max=0.0;
		for (j=0;j<n; j++){
			for (i=j+1; i<n; i++){
				if (max<fabs(A[i][j])){
					max = fabs(A[i][j]);
				}
			}
		}

	}//end loop:

	//imprimindo A:
	for (i=0; i<n; i++){
		for (j=0; j<n; j++){
			printf("A[%d][%d]=%f\t", i,j,A[i][j]);
		}
		printf("\n");
	}


	for (i=0; i<n; i++){
		free(U[i]);
		free(P[i]);
	}
	free(U);
	free(P);

}
