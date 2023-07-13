#include "GaussianElimination.h"


void GaussianElimination(double **A, double **S, double *b, int n){

	int k, i, j, p;
	double aux;

	k=0;
	//Guassian elimination - Triangularization
	for (i=0; i<n-1; i++){
		//k=0;
		if (A[i][i]==0.0){				
			printf("\n Need to change line \n");
			//Searching for a new first line
			while ((A[i+k][i]==0.0) && (k<n)){
				k++;
			}
			//changing lines of the matrix
			for (p=i; p<n; p++){
				aux = A[i][p];
				A[i][p] = A[i+k][p];
				A[i+k][p] = aux;
			}
			//change independent vector b.
			aux = b[i];
			b[i] = b[i+k];
			b[i+k] = aux;

		}

		k=i+1;			
		while ( k < n )
		{
			for (j=i; j<n; j++){
			
				//printf("A[k][j]= %f \n", A[k][j]);
				S[k][j] = A[k][j] - A[i][j]*(A[k][i]/A[i][i]);
			}
			b[k] = b[k] - (b[i]*(A[k][i]/A[i][i]));

			//Update results - isso me ajuda a evitar ter que trabalhar com a S na equação anterior.
			for (j=i; j<n; j++){
				A[k][j] = S[k][j];
				//printf("A[k][j]= %f \n", A[k][j]);
				//printf("here: k=%d, j=%d \n", k, j);getchar();

			}

			k++;
		}


	}//end for lines i



}
