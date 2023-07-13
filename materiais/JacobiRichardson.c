#include "JacobiRichardson.h"

void JacobiRichardson(double **A, double *s, double *b, int n){
	double B[n][n], g[n], s_old[n], sum, ItMax=1.0e+3;
	double tol=1.0e-14, error=1.0, aux1, aux2;
	int i,j, count=0;
	for(i=0;i<n;i++){
		s_old[i]=0.0;//inicialização de s
		g[i]=(1.0/A[i][i])*b[i];//g
		for(j=0;j<n;j++){
			if(i==j){
				B[i][j]=0.0;
			}else{
			B[i][j]= -(1.0/A[i][i])*A[i][j];
			}
		}
	}
	
	while(count<ItMax){
		count++;
	for(i=0;i<n;i++){
	sum=0.0;
		for(j=0;j<n;j++){			
			sum += B[i][j]*s_old[j];
		}
		s[i]= sum + g[i];
	}
	aux1=0.0;
	aux2=0.0;
	for(i=0;i<n;i++){
	aux1 += (s[i]-s_old[i])*(s[i]-s_old[i]);
	aux2 += s[i]*s[i];
	}

	error = sqrt(aux1)/sqrt(aux2);

	for(i=0;i<n;i++){
		s_old[i]=s[i];
	}

	}//end while



}
