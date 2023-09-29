#include "MetodoPotencia.h"

double MetodoPotencia(double **A, double *u, int n, double tol, int inv){

	//u eh o autovetor
	double lambda, lambdaOld, erro, erroOld;
	double *y, *yOld, *z, *zOld, **As, alpha, detA;
	int i, j, r, iter, iterMax=1000, lu;

	y = (double *)malloc(sizeof(double)*n);
	yOld = (double *)malloc(sizeof(double)*n);
	z = (double *)malloc(sizeof(double)*n);
	zOld = (double *)malloc(sizeof(double)*n);
	As = (double **)malloc(sizeof(double *)*n);
	for (i=0; i<n; i++){
		As[i] =(double *)malloc(sizeof(double *)*n);
	}

	//printf("\ninv=%d\n", inv);getchar();

	//valor inicial de y:
	for (i=0; i<n; i++){
		y[i] = 1.0;
		yOld[i] = 1.0;
	}

	if (inv == 0){
		z = MatrixProduct(A, y, n);
		zOld=MatrixProduct(A, yOld, n);
	}

	erro = 100.0;	
	iter = 0;
	while ( (erro > tol) && (iter < iterMax) ){

		for (i=0; i<n; i++){
			for (j=0; j<n; j++){
				As[i][j] = A[i][j];
			}
		}

		//copia valor novo para velho
		//erroOld = erro;
		lambdaOld = lambda;
		for (i=0; i<n; i++){
			yOld[i] = y[i];
			zOld[i] = z[i];
		}

		//Esolha do metodo da potencia ou potencia inversa
		if (inv == 0){
			z = MatrixProduct(A, y, n);			
		}else{
			lu=1;
			detA=LUdecomposition(As, n);
			SistemaTriangularInferior(As, y, z, n, lu);
			SistemaTriangularSuperior(As, z, z, n);//O vetor independente eh a solucao do sistema triangular inferior.
		}
		alpha = MaxNorm(z, n);		
		for (i=0; i<n; i++){
			y[i] = z[i]/alpha;
			//printf("z[%d]=%f\n", i, z[i]);//getchar();
			//printf("y[%d]=%f\n", i, y[i]);//getchar();
		}

		r=0;//reinicia indice que ira calcular lambda
		for (i=0; i<n; i++){
			if (inv == 0){
				lambda = z[i]/yOld[i];
				//printf("lambda=%f\n", lambda);//getchar();

			}else{
				lambda = fabs(z[i]/yOld[i]);
				//printf("lambda=%f\n", lambda);getchar();
			}

			//para evitar nan
			//if (lambda == 0.0){
			//	lambda = 1.0;
			//}
			erro = fabs(lambda - lambdaOld)/fabs(lambda);
			if (i==0){
				erroOld = erro;
			}
			
			//printf("lambda=%f\n", lambda);
			if (erro > erroOld){
				//printf("erro=%f, erroOld=%f\n", erro, erroOld);
				erro = erroOld;
				lambda = lambdaOld;				
			}else{
				r = i;
			}

			//printf("i=%d\n", i);
			//printf("y=%f, yOld=%f\n", y[i], yOld[i]);
			//printf("lambda=%f\n", lambda);
			//printf("erro=%f\n", erro);getchar();

		}
		iter++;
	}
	printf("lambda=%f\n", lambda);
	printf("iter=%d\n", iter);
	//atualizar autovetor
	for (i=0; i<n; i++){
		u[i] = y[i];
		free(As[i]);
	}

	free(y);
	free(yOld);
	free(z);
	free(zOld);
	free(As);


	if (inv == 0){
		return lambda;
	}else{
		return (1.0/lambda);
	}
}
