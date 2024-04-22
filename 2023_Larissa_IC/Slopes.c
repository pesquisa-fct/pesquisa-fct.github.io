#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main(){

	FILE *slopes;
	int N=100, i;
	double x0, xf, h, a, b;
	double *x, *u, *v, *p, *psi, *omega;

	//memory allocation
	x = (double *)malloc(sizeof(double)*(N+1));
	u = (double *)malloc(sizeof(double)*(N+1));
	v = (double *)malloc(sizeof(double)*(N+1));
	p = (double *)malloc(sizeof(double)*(N+1));
	psi = (double *)malloc(sizeof(double)*(N+1));
	omega = (double *)malloc(sizeof(double)*(N+1));

	//Parameters
	x0=-1.3;
	xf=-0.6;
	h = (xf-x0)/N;

	//Discretizacao do dominio:
	x[0] = x0;
	for (i=0; i<N; i++){
		x[i+1] = x[i] + h;
	}

	slopes = fopen("Slopes4Lx4L_M1_Re1.dat", "w");

	//Calculo dos slopes:
	for (i=0; i<=N; i++){		
		u[i] = (0.5445)*x[i] + 0.365;
		v[i] = (0.5445)*x[i] - 0.2;
		p[i] = (-0.4555)*x[i] + 2.1;
		psi[i] = (1.5445)*x[i] - 0.4;
		omega[i] = (-0.4555)*x[i] - 0.6;
		fprintf(slopes,"%f \t %f \t %f \t %f \t %f \t %f \n ", x[i], u[i], v[i], p[i], psi[i], omega[i]);
	}
	fclose(slopes);

	free(x);
	free(u);
	free(v);
	free(p);
	free(psi);
	free(omega);

	return (0);
}
