#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "GaussSeidel.h"
#include "ConjugateGradientMethod.h"

#include "FuncoesAuxiliares.h"
#include "Visualizacao.h"

typedef struct Indices {
	int i, j;
} indices;	

int main(){

	double t0=0.0, tf = 10.0, x0=0.0, xf=8.0, y0=0.0, yf=4.0;
	double xm=4.0, ymDown=1.5, ymUp=2.5;
	int Nt=100000, Nx=80, Ny=40;
	double *x, *y, *t; 
	double **uOld, **vOld, **psiOld, **omegaOld;
	double **uNew, **vNew, **psiNew, **omegaNew;
	double **pOld, **pNew;
	int **VecInd;
	double c, ConvU, ConvV, ViscX, ViscY;
	double dx, dy, dt;
	int i, j, n, k, p, dim;
	int i0, im, ie;		//i0: inicial; im: meio; ie: final
	int j0, jm1, jm2, je;
	indices *vec;
	double Re, beta2, aux;
	int LinearSystem;

	//Parameters:
	Re = 1.0;
	//_______________________________
	//Gauss-Seidel:			0
	//Conjugate Gradient Method:	1
	LinearSystem = 1;

	dx = (xf-x0)/Nx;
	dy = (yf-y0)/Ny;
	dt = (tf-t0)/Nt;

	//Memory allocation --------------------------------
	x = (double *)malloc(sizeof(double)*(Nx+1));
	y = (double *)malloc(sizeof(double)*(Ny+1));
	t = (double *)malloc(sizeof(double)*(Nt+1));
	//--------------------------------------------------
	uOld = (double **)malloc(sizeof(double *)*(Nx+1));
	vOld = (double **)malloc(sizeof(double *)*(Nx+1));
	psiOld = (double **)malloc(sizeof(double *)*(Nx+1));
	omegaOld = (double **)malloc(sizeof(double *)*(Nx+1));
	pOld = (double **)malloc(sizeof(double *)*(Nx+1));
	//--------------------------------------------------
	uNew = (double **)malloc(sizeof(double *)*(Nx+1));
	vNew = (double **)malloc(sizeof(double *)*(Nx+1));
	psiNew = (double **)malloc(sizeof(double *)*(Nx+1));
	omegaNew = (double **)malloc(sizeof(double *)*(Nx+1));
	pNew = (double **)malloc(sizeof(double *)*(Nx+1));

	VecInd = (int **)malloc(sizeof(int *)*(Nx + 1));
	for (i=0; i<=Nx; i++){
		uOld[i] = (double *)malloc(sizeof(double)*(Ny+1));
		vOld[i] = (double *)malloc(sizeof(double)*(Ny+1));
		psiOld[i] = (double *)malloc(sizeof(double)*(Ny+1));
		omegaOld[i] = (double *)malloc(sizeof(double)*(Ny+1));
		pOld[i] = (double *)malloc(sizeof(double)*(Ny+1));
		//--------------------------------------------------
		uNew[i] = (double *)malloc(sizeof(double)*(Ny+1));
		vNew[i] = (double *)malloc(sizeof(double)*(Ny+1));
		psiNew[i] = (double *)malloc(sizeof(double)*(Ny+1));
		omegaNew[i] = (double *)malloc(sizeof(double)*(Ny+1));
		pNew[i] = (double *)malloc(sizeof(double)*(Ny+1));

		VecInd[i] = (int *)malloc(sizeof(int)*(Ny + 1));
	}

	//Discretizacao do dominio:
	x[0] = x0;
	for (i=0; i<Nx; i++){
		x[i+1] = x[i] + dx;
	}
	y[0] = y0;
	for (j=0; j<Ny; j++){
		y[j+1] = y[j] + dy;
	}
	t[0] = t0;
	for (n=0; n<Nt; n++){
		t[n+1] = t[n] + dt;
	}
	printf("\n Apos discretizacao do dominio \n");

	//-------------------------------------
	n=0;
	//Inicializacao:
	for (i=0; i<=Nx; i++){
		for (j=0; j<=Ny; j++){
			uOld[i][j] = uNew[i][j] = 0.0;
			vOld[i][j] = vNew[i][j] = 0.0;
			pOld[i][j] = pNew[i][j] = 0.0;
			psiOld[i][j] = psiNew[i][j] = 0.0;
			omegaOld[i][j] = omegaNew[i][j] = 0.0;
		}
	}
	printf("\n apos inicializacao Ny=%d\n\n \n", Ny);
	//Boundary conditions:
	//INFLOW:
	i=0;
	//constntant c:
	c = 0.0;
	for (j=0; j<=Ny; j++){
		//printf("j=%d\n", j);
		uOld[i][j] = uNew[i][j] = -(3.0/32.0)*( y[j] - 4.0 )*y[j];
		vOld[i][j] = vNew[i][j] = 0.0;
		psiOld[i][j] = psiNew[i][j] = -(1.0/32.0)*y[j]*y[j]*y[j] + (3.0/16.0)*y[j]*y[j] + c;
		//psiOld[i][j] = psiOld[i+1][j] + dy*uOld[i][j];
		//psiNew[i][j] = psiNew[i+1][j] + dy*uNew[i][j];
		omegaOld[i][j] = omegaNew[i][j] = (3.0/16.0)*y[j] - (3.0/8.0);
	}
	//PAREDES:
	
	//OUTFLOW:
	//-------------------------------------
	printf("\n Apos condicoes de contorno \n");
	for (i=0; i<=Nx; i++){
		//x0:
		if ( fabs(x[i] - x0) < 1.0e-06 ){
			i0 = i;
		}
		//xm:
		if ( fabs(x[i] - xm) < 1.0e-06 ){
			im = i;
		}
		//xf:
		if ( fabs(x[i] - xf) < 1.0e-06 ){
			ie = i;
		}
	}
	//printf("i0=%d, \t im=%d \t ie=%d\n", i0, im, ie);
	//printf("x[i0]=%f, \t x[im]=%f, \t x[ie]=%f\n", x[i0], x[im], x[ie]);getchar();

	for (j=0; j<=Ny; j++){
		//y0:
		if ( fabs(y[j] - y0) < 1.0e-06 ){
			j0 = j;
		}
		//ym1:
		if ( fabs(y[j] - ymDown) < 1.0e-06 ){
			jm1 = j;
		}
		//ym2:
		if ( fabs(y[j] - ymUp) < 1.0e-06 ){
			jm2 = j;
		}
		//ye:
		if ( fabs(y[j] - yf) < 1.0e-06 ){
			je = j;
		}
	}
	//printf("y[j0]=%f, \t y[jm1]=%f, \t y[jm2]=%f, \t y[je]=%f \n", y[j0], y[jm1], y[jm2], y[je]);getchar();

	//_______________________________________________________________________
	//Alocando memoria para os indices:
	//_______________________________________________________________________
	p=0;
	for (i=i0+1; i<im; i++){
		for (j=j0+1; j<je; j++){
			p++;
		}
	}
	//bloco 2:
	for (i=im; i<ie; i++){
		for (j=jm1+1; j<jm2; j++){
			p++;
		}
	}
	p=p-1;
	//_______________________________________________________________________
	vec = (indices *)malloc(sizeof(indices)*(p+1));

	//________________________________________________________________________
	//Armazenando indices i,j e relacionando com a posicao do vetor solucao p:
	//bloco 1:
	p=0;
	for (i=i0+1; i<im; i++){
		for (j=j0+1; j<je; j++){
			vec[p].i = i;
			vec[p].j = j;
			VecInd[i][j] = p;
			p++;
		}
	}
	//printf("vec[p-1].i=%d, vec[p-1].j=%d\n", vec[p-1].i, vec[p-1].j);
	//printf("p=%d\n", p);getchar();
	//bloco 2:
	for (i=im; i<ie; i++){
		for (j=jm1+1; j<jm2; j++){
			vec[p].i = i;
			vec[p].j = j;
			VecInd[i][j] = p;
			p++;
		}
	}
	dim = p;
	p=p-1;
	//_______________________________________________________________________
	//printf("dim=%d\n", dim);getchar();
	
	/*
	for (k=0; k<dim; k++){
		printf("k=%d\n", k);
		printf("i=%d, j=%d\n", vec[k].i, vec[k].j);//getchar();
	}
	*/
	


	//Dados para a resolução da equação de Poisson:
	double **A, *s, *b;
	int iter;
	A=(double **)malloc(sizeof(double *)*dim);
	for (k=0; k<dim; k++){
		A[k]=(double *)malloc(sizeof(double)*dim);
	}
	b = (double *)malloc(sizeof(double)*dim);
	s = (double *)malloc(sizeof(double)*dim);

	//Salvando estado inicial
	//if (n % 100 == 0.0){
	ImprimeArquivoVTK(uNew, vNew, pNew, psiNew, omegaNew, x0, xf, y0, yf, x, y, Nx, Ny, n);
	//}

	while (n*dt < tf){
		n++;

		//____________________________________________________________________________
		//Atualiza condicoes de contorno para vorticidade:
		//printf("Atualiza condicoes de contorno:\n");
		//_______
		//INFLOW:
		i=i0;
		for (j=0; j<=je; j++){
			//uOld[i][j] = uNew[i][j] = -(3.0/32.0)*( y[j] - 4.0 )*y[j];
			//vOld[i][j] = vNew[i][j] = 0.0;
			//psiOld[i][j] = psiOld[i+1][j] + dy*uOld[i][j];
			//omegaOld[i][j] = (-85.0*psiOld[i][j] + 108.0*psiOld[i+1][j] - 27.0*psiOld[i+2][j] + 4.0*psiOld[i+3][j])/(18.0*dx*dx);
			//psiOld[i][j] = psiOld[i+1][j] + dy*uOld[i][j];
			//psiNew[i][j] = psiNew[i+1][j] + dy*uOld[i][j];

		}
		//printf("INFLOW\n");
		//_______________
		//PAREDE DIREITA:
		i=im;
		//printf("omegaOld[][]=%f\n", omegaOld[1][1]);
		for (j=0; j<=jm1; j++){
			//printf("im=%d, j=%d\n", i, j);
			omegaNew[i][j] = omegaOld[i][j] = (-85.0*psiOld[i][j] + 108.0*psiOld[i-1][j] - 27.0*psiOld[i-2][j] + 4.0*psiOld[i-3][j])/(18.0*dx*dx);
			//omegaOld[i][j] = -(psiOld[i][j] - 2.0*psiOld[i-1][j] + psiOld[i-2][j])/(dx*dx);
			//printf("omegaOld[%d][%d]=%f\n", i, j, omegaOld[i][j]);
		}
		for (j=jm2; j<=je; j++){
			omegaNew[i][j] = omegaOld[i][j] = (-85.0*psiOld[i][j] + 108.0*psiOld[i-1][j] - 27.0*psiOld[i-2][j] + 4.0*psiOld[i-3][j])/(18.0*dx*dx);
			//omegaOld[i][j] = -(psiOld[i][j] - 2.0*psiOld[i-1][j] + psiOld[i-2][j])/(dx*dx);
		}
		//printf("PAREDE DIREITA\n");
		//_____________
		//PAREDE BAIXO:
		j=j0;
		for (i=0; i<=im; i++){
			omegaNew[i][j] = omegaOld[i][j] = (-85.0*psiOld[i][j] + 108.0*psiOld[i][j+1] - 27.0*psiOld[i][j+2] + 4.0*psiOld[i][j+3])/(18.0*dy*dy);
			//omegaOld[i][j] = - (psiOld[i][j] - 2.0*psiOld[i][j+1] + psiOld[i][j+2])/(dy*dy);
			omegaNew[i][j] = omegaOld[i][j] = (1.0/( 18.0*(dy*dy)) )*( -85.0*psiOld[i][j] + 108.0*psiOld[i][j+1] - 27.0*psiOld[i][j+2] + 4.0*psiOld[i][j+3] );//Goiabinha
		}
		j=jm1;
		for (i=im; i<=ie; i++){
			omegaOld[i][j] = (-85.0*psiOld[i][j] + 108.0*psiOld[i][j+1] - 27.0*psiOld[i][j+2] + 4.0*psiOld[i][j+3])/(18.0*dy*dy);
			//omegaOld[i][j] = - (psiOld[i][j] - 2.0*psiOld[i][j+1] + psiOld[i][j+2])/(dy*dy);
		}
		//printf("PAREDE BAIXO\n");
		//________________
		//PAREDE SUPERIOR:
		j=je;
		for (i=0; i<=im; i++){
			omegaNew[i][j] = omegaOld[i][j] = (-85.0*psiOld[i][j] + 108.0*psiOld[i][j-1] - 27.0*psiOld[i][j-2] + 4.0*psiOld[i][j-3])/(18.0*dy*dy);
			//omegaOld[i][j] = - (psiOld[i][j] - 2.0*psiOld[i][j-1] + psiOld[i][j-2])/(dy*dy);
		}
		j=jm2;
		for (i=im; i<=ie; i++){
			omegaNew[i][j] = omegaOld[i][j] = (-85.0*psiOld[i][j] + 108.0*psiOld[i][j-1] - 27.0*psiOld[i][j-2] + 4.0*psiOld[i][j-3])/(18.0*dy*dy);
			//omegaOld[i][j] = - (psiOld[i][j] - 2.0*psiOld[i][j-1] + psiOld[i][j-2])/(dy*dy);
		}
		//printf("PAREDE SUPERIOR\n");
		//________
		//OUTFLOW:
		i=ie;
		for (j=jm1; j<=jm2; j++){
			omegaOld[i][j] = omegaOld[i-1][j];
			psiOld[i][j] = psiOld[i-1][j];
		}
		//printf("OUTFLOW\n");
		//----------------------------------------------------------------------------

		//____________________________________________________________________________
		//bloco 1:
		//printf("bloco 1:\n");
		for (i=1; i<im; i++){
			for (j=1; j<je; j++){

				//INFLOW:
				if (i==1){
					//ConvU = (omegaOld[i+1][j] - omegaOld[i][j])/dx;
					ConvU = -uOld[i][j]*(omegaOld[i+1][j] - omegaOld[i-1][j])/(2.0*dx);//valor de omega eh conhecido no inflow
					/*									   
					//UpWind													   
					if (uOld[i][j] >= 0.0){
						ConvU = -uOld[i][j]*(omegaOld[i][j] - omegaOld[i-1][j])/(dx);
					}else{
						ConvU = -uOld[i][j]*(omegaOld[i+1][j] - omegaOld[i][j])/(dx);
					}
					*/								   
					ViscX = (omegaOld[i+1][j] - 2.0*omegaOld[i][j] + omegaOld[i-1][j])/(dx*dx);//Posso usar diferenca centrada.
				//PAREDE DIREITA:											       
				}else if ( (i==im-1) && ((j<=jm1) || (j>=jm2)) ){					
					//ConvU = -u[i][j](omegaOld[i][j] - omegaOld[i-1][j])/(dx);//Diferenca regressiva
					ConvU = -uOld[i][j]*(omegaOld[i+1][j] - omegaOld[i-1][j])/(2.0*dx);//Diferença centrada			 
					/*													   
					//UpWind													   
					if (uOld[i][j] >= 0.0){
						ConvU = -uOld[i][j]*(omegaOld[i][j] - omegaOld[i-1][j])/(dx);
					}else{
						ConvU = -uOld[i][j]*(omegaOld[i+1][j] - omegaOld[i][j])/(dx);
					}
					*/
					ViscX = (omegaOld[i+1][j] - 2.0*omegaOld[i][j] + omegaOld[i-1][j])/(dx*dx);//Diferenca centrada, pois omegaOld[i+1][j] conheco.
				}else{
					ConvU = -uOld[i][j]*(omegaOld[i+1][j] - omegaOld[i-1][j])/(2.0*dx);//Diferenca centrada					
					/*
					//UpWind													   
					if (uOld[i][j] >= 0.0){
						ConvU = -uOld[i][j]*(omegaOld[i][j] - omegaOld[i-1][j])/(dx);
					}else{
						ConvU = -uOld[i][j]*(omegaOld[i+1][j] - omegaOld[i][j])/(dx);
					}
					*/
					ViscX = (omegaOld[i+1][j] - 2.0*omegaOld[i][j] + omegaOld[i-1][j])/(dx*dx);
				}
				//PAREDE CIMA:
				if (j==je-1){
					ConvV = -vOld[i][j]*(omegaOld[i][j+1] - omegaOld[i][j-1])/(2.0*dy);
					/*
					//UpWind													   
					if (vOld[i][j] >= 0.0){
						ConvV = -vOld[i][j]*(omegaOld[i][j] - omegaOld[i][j-1])/(dy);
					}else{
						ConvV = -vOld[i][j]*(omegaOld[i][j+1] - omegaOld[i][j])/(dy);
					}
					*/
					ViscY = (omegaOld[i][j+1] - 2.0*omegaOld[i][j] + omegaOld[i][j-1])/(dy*dy);//Posso usar diferenca centrada.
				//PAREDE BAIXO:
				}else if (j==1){
					ConvV = -vOld[i][j]*(omegaOld[i][j+1] - omegaOld[i][j-1])/(2.0*dy);//Diferenca progressiva
					/*
					//UpWind													   
					if (vOld[i][j] >= 0.0){
						ConvV = -vOld[i][j]*(omegaOld[i][j] - omegaOld[i][j-1])/(dy);
					}else{
						ConvV = -vOld[i][j]*(omegaOld[i][j+1] - omegaOld[i][j])/(dy);
					}
					*/
					ViscY = (omegaOld[i][j+1] - 2.0*omegaOld[i][j] + omegaOld[i][j-1])/(dy*dy);//Posso usar diferenca centrada.
				}else{
					ConvV = -vOld[i][j]*(omegaOld[i][j+1] - omegaOld[i][j-1])/(2.0*dy);//Diferenca centrada
					/*
					//UpWind													   
					if (vOld[i][j] >= 0.0){
						ConvV = -vOld[i][j]*(omegaOld[i][j] - omegaOld[i][j-1])/(dy);
					}else{
						ConvV = -vOld[i][j]*(omegaOld[i][j+1] - omegaOld[i][j])/(dy);
					}
					*/
					ViscY = (omegaOld[i][j+1] - 2.0*omegaOld[i][j] + omegaOld[i][j-1])/(dy*dy);//Posso usar diferenca centrada.
				}

				//Euler Explicito - omegaNew
				//omegaNew[i][j] = omegaOld[i][j] + dt*( ConvU + ConvV + (1.0/Re)*( ViscX + ViscY ) );
				omegaNew[i][j] = omegaOld[i][j] + dt*( -uOld[i][j]*(omegaOld[i+1][j] - omegaOld[i-1][j])/(2.0*dx) - vOld[i][j]*(omegaOld[i][j+1] - omegaOld[i][j-1])/(2.0*dy) + (1.0/Re)*( (omegaOld[i+1][j] - 2.0*omegaOld[i][j] + omegaOld[i-1][j])/(dx*dx) + (omegaOld[i][j+1] - 2.0*omegaOld[i][j] + omegaOld[i][j-1])/(dy*dy) ) );
			}
		}
		//bloco 2:
		for (i=im; i<ie; i++){
			for (j=jm1+1; j<jm2; j++){
				
				//OUTFLOW
				if (i==ie-1){
					ConvU = -uOld[i][j]*(omegaOld[i+1][j] - omegaOld[i-1][j])/(2.0*dx);
					/*
					//UpWind													   
					if (uOld[i][j] >= 0.0){
						ConvU = -uOld[i][j]*(omegaOld[i][j] - omegaOld[i-1][j])/(dx);
					}else{
						ConvU = -uOld[i][j]*(omegaOld[i+1][j] - omegaOld[i][j])/(dx);
					}
					*/
					ViscX = (omegaOld[i+1][j] - 2.0*omegaOld[i][j] + omegaOld[i-1][j])/(dx*dx);//Diferenca centrada, pois omegaOld[i+1][j] conheco.
				}else{
					ConvU = -uOld[i][j]*(omegaOld[i+1][j] - omegaOld[i-1][j])/(2.0*dx);
					/*
					//UpWind													   
					if (uOld[i][j] >= 0.0){
						ConvU = -uOld[i][j]*(omegaOld[i][j] - omegaOld[i-1][j])/(dx);
					}else{
						ConvU = -uOld[i][j]*(omegaOld[i+1][j] - omegaOld[i][j])/(dx);
					}
					*/
					ViscX = (omegaOld[i+1][j] - 2.0*omegaOld[i][j] + omegaOld[i-1][j])/(dx*dx);//Diferenca centrada, pois omegaOld[i+1][j] conheco.
				}
				//PAREDE CIMA:
				if (j==jm2-1){
					//ConvV = -v[i][j](omegaOld[i][j] - omegaOld[i][j-1])/(dy);//Diferenca regressiva
					ConvV = -vOld[i][j]*(omegaOld[i][j+1] - omegaOld[i][j-1])/(2.0*dy);//Diferenca centrada
					/*
					//UpWind													   
					if (vOld[i][j] >= 0.0){
						ConvV = -vOld[i][j]*(omegaOld[i][j] - omegaOld[i][j-1])/(dy);
					}else{
						ConvV = -vOld[i][j]*(omegaOld[i][j+1] - omegaOld[i][j])/(dy);
					}
					*/
					ViscY = (omegaOld[i][j+1] - 2.0*omegaOld[i][j] + omegaOld[i][j-1])/(dy*dy);//Posso usar diferenca centrada.
				//PAREDE BAIXO:
				}else if (j==jm1+1){
					//ConvV = -vOld[i][j](omegaOld[i][j+1] - omegaOld[i][j])/(dy);//Diferenca progressiva
					ConvV = -vOld[i][j]*(omegaOld[i][j+1] - omegaOld[i][j-1])/(2.0*dy);//Diferenca centrada
					/*
					//UpWind													   
					if (vOld[i][j] >= 0.0){
						ConvV = -vOld[i][j]*(omegaOld[i][j] - omegaOld[i][j-1])/(dy);
					}else{
						ConvV = -vOld[i][j]*(omegaOld[i][j+1] - omegaOld[i][j])/(dy);
					}
					*/
					ViscY = (omegaOld[i][j+1] - 2.0*omegaOld[i][j] + omegaOld[i][j-1])/(dy*dy);//Posso usar diferenca centrada.
				}else{
					ConvV = -vOld[i][j]*(omegaOld[i][j+1] - omegaOld[i][j-1])/(2.0*dy);//Diferenca centrada
					/*
					//UpWind													   
					if (vOld[i][j] >= 0.0){
						ConvV = -vOld[i][j]*(omegaOld[i][j] - omegaOld[i][j-1])/(dy);
					}else{
						ConvV = -vOld[i][j]*(omegaOld[i][j+1] - omegaOld[i][j])/(dy);
					}
					*/
					ViscY = (omegaOld[i][j+1] - 2.0*omegaOld[i][j] + omegaOld[i][j-1])/(dy*dy);//Posso usar diferenca centrada.
				}

				//Euler Explicito - omegaNew
				//omegaNew[i][j] = omegaOld[i][j] + dt*( ConvU + ConvV + (1.0/Re)*( ViscX + ViscY ) );
				omegaNew[i][j] = omegaOld[i][j] + dt*( -uOld[i][j]*(omegaOld[i+1][j] - omegaOld[i-1][j])/(2.0*dx) - vOld[i][j]*(omegaOld[i][j+1] - omegaOld[i][j-1])/(2.0*dy) + (1.0/Re)*( (omegaOld[i+1][j] - 2.0*omegaOld[i][j] + omegaOld[i-1][j])/(dx*dx) + (omegaOld[i][j+1] - 2.0*omegaOld[i][j] + omegaOld[i][j-1])/(dy*dy) ) );

			}
		}


		//################################################################
		//Sol. Eq. Poisson:
		//################################################################
		//_______________________________________
		//Atualiza condição de contorno para psi:		
		//PAREDE VERTICAL ESQUERDA (INFLOW):
		i=i0;
		for (j=0; j<=je; j++){
			//psiOld[i][j] = 2.0*psiOld[i+1][j]  - psiOld[i+2][j];
		}
		
		//PAREDE HORIZONTAL SUPERIOR:
		j=je;
		for (i=1; i<=im; i++){
			psiOld[i][j] = psiNew[i][j] = y[j]/4.0;
			//psiOld[i][j] = 2.0*psiOld[i][j-1] - 1.0*psiOld[i][j-2];
			//psiNew[i][j] = 2.0*psiNew[i][j-1] - 1.0*psiNew[i][j-2];

		}
		j=jm2;
		for (i=im; i<=ie; i++){
			psiOld[i][j] = psiNew[i][j] = 1.0;//y[j]/4.0;
			//psiOld[i][j] = psiOld[i][j-1];
			//psiNew[i][j] = psiNew[i][j-1];

		}
		//PAREDE HORIZONTAL INFERIOR:
		j=j0;
		for (i=1; i<=im; i++){
			psiOld[i][j] = psiNew[i][j] = y[j]/4.0;
			//psiOld[i][j] = psiOld[i][j+1];
			//psiNew[i][j] = psiNew[i][j+1];
		}
		j=jm1;
		for (i=im; i<=ie; i++){
			psiOld[i][j] = psiNew[i][j] = 0.0;//y[j]/4.0;
			//psiOld[i][j] = psiOld[i][j+1];
			//psiNew[i][j] = psiNew[i][j+1];
		}
		//PAREDE VERTICAL DIREITA:
		i=im;
		for (j=0; j<=jm1; j++){
			psiOld[i][j] = psiNew[i][j] = 0.0;//y[j]/4.0;
			//psiOld[i][j] = psiOld[i-1][j];
			//psiNew[i][j] = psiNew[i-1][j];
		}
		for (j=jm2; j<=je; j++){
			psiOld[i][j] = psiNew[i][j] = 1.0;//y[j]/4.0;
			//psiOld[i][j] = psiOld[i-1][j];
			//psiNew[i][j] = psiNew[i-1][j];
		}
		//PAREDE VERTICAL DIREITA - OUTFLOW: Neumann Homogeneo
		i=ie;
		for (j=jm1; j<=jm2; j++){
			psiOld[i][j] = psiOld[i-1][j];
			psiNew[i][j] = psiNew[i-1][j];
		}
		//________________________________________________________________
		//Inicializa matriz dos coeficientes A, vetor b e chute inicial (tbm solução) s:
		for (i=0; i<dim; i++){
			b[i] = 0.0;
			s[i] = 0.0;
			for (j=0; j<dim; j++){
				A[i][j] = 0.0;
			}
		}

		/*		
		//_______________________________________________________
		//Monta matriz dos coeficientes A e vetor independente b:
		beta2 = (dx/dy)*(dx/dy);
		aux = 1.0/( 2.0*(1 + beta2) );
		for (p=0; p<dim; p++){
			//indices:
			i = vec[p].i;
			j = vec[p].j;

			//printf("x[%d]=%f, y[%d]=%f\n", i, x[i], j, y[j]);getchar();

			//Diagonal principal e vetor b:
			//A[p][p] = - (2.0/(dx*dx)) - (2.0/(dy*dy));
			A[p][p] = 1.0;
			b[p] = -aux*dx*dx*omegaNew[i][j];

			//PAREDE VERTICAL ESQUERDA (INFLOW):
			if ( fabs(x[i-1] - x0)<1.0e-06){
				//printf("parede lateral esquerda\n");
				//printf("i=%d, j=%d\n", i, j);
				//printf("x=%f, y=%f\n", x[i], y[j]);getchar();

				//b[p] += - (psiOld[i-1][j]/(dx*dx));
				b[p] += aux*psiOld[i-1][j];
			}else{
				k = VecInd[i-1][j];
				//A[p][k] = 1.0/(dx*dx);
				A[p][k] = -aux;
			}
			//PAREDE HORIZONTAL SUPERIOR (SOLID WALL):
			if ( ( fabs(y[j+1] - yf)<1.0e-06 && (x[i]>x0 && x[i]<xm) ) || ( fabs(y[j+1] - ymUp)<1.0e-06 && ( (x[i]>xm || fabs(x[i]-xm)<1.0e-06) && x[i]<xf)) ){
				//printf("parede superior\n");
				//printf("i=%d, j=%d\n", i, j);
				//printf("x=%f, y=%f\n", x[i], y[j]);getchar();
				//b[p] += - (psiOld[i][j+1]/(dy*dy));
				b[p] += beta2*aux*psiOld[i][j+1];
			}else{
				//printf("y[j+1]=%f, yf=%f\n", y[j+1], yf);
				//printf("i=%d, j=%d\n", i, j);
				//printf("x=%f, y=%f\n", x[i], y[j]);getchar();
				//printf("aqui\n");getchar();
				//A[p][p+1] = 1.0/(dy*dy);
				A[p][p+1] = -beta2*aux;
			}
			//PAREDE HORIZONTAL INFERIOR (SOLID WALL):
			if ( ( fabs(y[j-1] - y0)<1.0e-06 && (x[i]>x0 && x[i]<xm) ) || ( fabs(y[j-1] - ymDown)<1.0e-06 && ( (x[i]>xm || fabs(x[i]-xm)<1.0e-06) && x[i]<xf)) ){
				//printf("parede inferior\n");
				//printf("i=%d, j=%d\n", i, j);
				//printf("x=%f, y=%f\n", x[i], y[j]);getchar();
				//b[p] += - (psiOld[i][j-1]/(dy*dy));
				b[p] += beta2*aux*psiOld[i][j-1];
			}else{
				//A[p][p-1] = 1.0/(dy*dy);
				A[p][p-1] = -beta2*aux;
			}
			//PAREDE LATERAL DIREITA (SOLID WALL) || OUTFLOW:
			if ( ( ( fabs(x[i+1] - xm)<1.0e-06 && ((y[j] > ymUp || fabs(y[j] - ymUp)<1.0e-06) || (y[j] < ymDown || fabs(y[j]-ymDown)<1.0e-06)) )  ) || (fabs(x[i+1]-xf)<1.0e-06) ){
				//printf("parede lateral direita\n");
				//printf("i=%d, j=%d\n", i, j);
				//printf("x=%f, y=%f\n", x[i], y[j]);getchar();

				//b[p] += - (psiOld[i+1][j]/(dx*dx));
				b[p] += aux*psiOld[i+1][j];
			}else{
				k = VecInd[i+1][j];
				//A[p][k] = 1.0/(dx*dx);
				A[p][k] = -aux; 
			}
		}
		*/
		/*
		//________________
		//Imprime matriz A:
		for (i=0; i<dim; i++){
			b[i] = 0.0;
			for (j=0; j<dim; j++){
				printf("A[%d][%d]=%f", i, j, A[i][j]);
			}
			printf("\n");getchar();
		}
		*/
		/*
		//_______________________________
		//Solução do sistema linear:
		//_______________________________
		switch (LinearSystem){
			case 0:
				//Gauss-Seidel:
				iter = GaussSeidel(A, b, s, dim);
			break;
			case 1:
				//Conjugate Gradient Method
				iter = ConjugateGradientMethod(A, b, s, dim, 1.0e-06);
			break;
			default:
				printf("Escolha um metodo válido para a resolução do sistema linear\n");
				exit(1);
			break;
		}
		*/

			
		//Explicit calculation
		beta2 = (dx/dy)*(dx/dy);
		aux = 1.0/( 2.0*(1 + beta2) );
		for (p=0; p<dim; p++){
			//printf("%f\n", s[p]);getchar();
			i=vec[p].i;
			j=vec[p].j;

			psiNew[i][j] = aux * ( psiOld[i - 1][j] + psiOld[i + 1][j] + (beta2) * ( psiOld[i][j + 1] + psiOld[i][j - 1] ) - (dx*dx) * omegaNew[i][j]);
		}
		
		//_______________________________________
		//Copia vetor s para a matriz psiNew:
		for (p=0; p<dim; p++){
			//printf("%f\n", s[p]);getchar();
			i=vec[p].i;
			j=vec[p].j;
			//psiNew[i][j] = s[p];

			//Calculo das velocidades:
			//uNew[i][j] = (psiNew[i][j+1] - psiNew[i][j-1])/(2.0*dy);
			//vNew[i][j] = -(psiNew[i+1][j] - psiNew[i-1][j])/(2.0*dx);
			
		  	uNew[i][j] =   (psiNew[i][j+1] - psiNew[i][j-1])/(2.0*dy);
                        vNew[i][j] = - (psiNew[i+1][j] - psiNew[i-1][j])/(2.0*dx);

			//Copia solução nova para velha:
			omegaOld[i][j] = omegaNew[i][j];
			psiOld[i][j] = psiNew[i][j];
			uOld[i][j] = uNew[i][j];
			vOld[i][j] = vNew[i][j];	
			pOld[i][j] = pNew[i][j];

		}
		/*
		for (i=1; i<ie; i++){
			for (j=1; j<je; j++){
				//Copia solução nova para velha:
				omegaOld[i][j] = omegaNew[i][j];
				psiOld[i][j] = psiNew[i][j];
				uOld[i][j] = uNew[i][j];
				vOld[i][j] = vNew[i][j];	
				pOld[i][j] = pNew[i][j];
			}
		}
		*/
		printf("n=%d\n", n);
		//Salvando valores:
		if (n % 10000 == 0.0){
			ImprimeArquivoVTK(uNew, vNew, pNew, psiNew, omegaNew, x0, xf, y0, yf, x, y, Nx, Ny, n);
		}

	}//end while

	//Vetores da equação de Poisson:
	for (k=0; k<dim; k++)	{
		free(A[k]);
	}
	free(A);
	free(b);
	free(s);

	//free memory
	free(vec);
	free(x);
	free(y);
	free(t);
	for (i=0; i<=Nx; i++){
		free(uOld[i]);
		free(vOld[i]);
		free(psiOld[i]);
		free(omegaOld[i]);
		free(pOld[i]);

		free(uNew[i]);
		free(vNew[i]);
		free(psiNew[i]);
		free(omegaNew[i]);
		free(pNew[i]);

		free(VecInd[i]);
	}
	free(uOld);
	free(vOld);
	free(psiOld);
	free(omegaOld);
	free(pOld);
	free(uNew);
	free(vNew);
	free(psiNew);
	free(omegaNew);		
	free(pNew);

	free(VecInd);

	vec=NULL;
	x=NULL;
	y=NULL;
	t=NULL;
	uOld=NULL;
	vOld=NULL;
	psiOld=NULL;
	omegaOld=NULL;
	pOld=NULL;

	uNew=NULL;
	vNew=NULL;
	psiNew=NULL;
	omegaNew=NULL;
	pNew=NULL;

	VecInd=NULL;

	return (0);
}