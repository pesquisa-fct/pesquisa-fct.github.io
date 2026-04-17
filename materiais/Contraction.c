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

	double t0, tf, x0, xf, y0, yf;
	double xm, ymDown, ymUp, ym, aInlet, Umax;
	int Nt, Nx, Ny;
	int Geometria, Malha;

	//-------------------------------------//						yf
	//					/
	//					/
	//					/
	//					/______________________________________//	ymUp
	//										/
	//										/
	//										/
	//					/---------------------------------------/	ymDown
	//					/
	//					/
	//					/
	//-------------------------------------//						y0										

	//x0					xm					xf

	//Informacoes do tempo:
	t0 = 0.0;
	tf = 10.0;
	Nt = 100000;
	//Dimensao da contracao 4:1
	//4L x 4L:	0
	//8L x 8L:	1
	//16L x 16L:	2
	//Dimensao da contracao 3:1
	//8L x 8L:	3
	//Dimensao da contracao 5:1
	//8L x 8L:	4
	Geometria = 1;
	//Malhas:
	//M1 (dx=dy=0.1):	0
	//M2 (dx=dy=0.05):	1
	//M3 (dx=dy=0.025):	2
	//M4 (dx=dy=0.0125):	3
	Malha = 1;
	switch (Geometria){
		case 0:
			//4L x 4L
			x0=0.0;
			xm=4.0;
			xf=8.0; 
			y0=0.0;
		        //ym=2.0;	
			ymDown=1.5;
			ymUp=2.5;
			yf=4.0;
			ym = 0.5*(yf+y0);
			aInlet = -3.0/8.0;
			switch(Malha){
				case 0:
					//Malha M1: dx=dy=0.1
					Nx=80; 
					Ny=40;
				break;
				case 1:
					//Malha M2: dx=dy=0.05
					Nx=160;
					Ny=80;
				break;
				case 2:
					//Malha M3: dx=dy=0.025
					Nx=320;
					Ny=160;
				break;
				case 3:
					//Malha M4: dx=dy=0.0125
					Nx=640;
					Ny=320;
				break;
				default:
					printf("Selecionar a opcao de malha corretamente\n");
					exit(1);
				break;
			}
		break;
		case 1:
			//8L x 8L
			x0=0.0;
			xm=8.0;
			xf=16.0; 
			y0=0.0;
		        ym=4.0;	
			ymDown=3.0;
			ymUp=5.0;
			yf=8.0;
			ym = 0.5*(yf+y0);
			aInlet = -3.0/8.0;
			Umax = 3.0/8.0;
			switch(Malha){
				case 0:
					//Malha M1: dx=dy=0.1
					Nx=160; 
					Ny=80;
				break;
				case 1:
					//Malha M2: dx=dy=0.05
					Nx=320;
					Ny=160;
				break;
				case 2:
					//Malha M3: dx=dy=0.025
					Nx=640;
					Ny=320;
				break;
				case 3:
					//Malha M4: dx=dy=0.0125
					Nx=1280;
					Ny=640;
				break;
				default:
					printf("Selecionar a opcao de malha corretamente\n");
					exit(1);
				break;
			}
		break;
		case 2:
			//16L x 16L
			x0=0.0;
			xm=16.0;
			xf=32.0; 
			y0=0.0; 
			//ym=4.0;
			ymDown=3.0;
			ymUp=5.0;
			yf=8.0;
			ym = 0.5*(yf+y0);
			aInlet = -3.0/8.0;
			switch(Malha){
				case 0:
					//Malha M1: dx=dy=0.1
					Nx=320; 
					Ny=80;
				break;
				case 1:
					//Malha M2: dx=dy=0.05
					Nx=640;
					Ny=160;
				break;
				case 2:
					//Malha M3: dx=dy=0.025
					Nx=1280;
					Ny=320;
				break;
				case 3:
					//Malha M4: dx=dy=0.0125
					Nx=2560;
					Ny=640;
				break;
				default:
					printf("Selecionar a opcao de malha corretamente\n");
					exit(1);
				break;
			}

		break;
		case 3://Contracao 3:1
			//8L x 8L
			x0=0.0;
			xm=8.0;
			xf=16.0; 
			y0=0.0;
		        //ym=4.0;	
			ymDown=2.0;
			ymUp=4.0;
			yf=6.0;
			ym = 0.5*(yf+y0);
			aInlet = -1.0/2.0;
			switch(Malha){
				case 0:
					//Malha M1: dx=dy=0.1
					Nx=160; 
					Ny=60;
				break;
				case 1:
					//Malha M2: dx=dy=0.05
					Nx=320;
					Ny=120;
				break;
				case 2:
					//Malha M3: dx=dy=0.025
					Nx=640;
					Ny=240;
				break;
				case 3:
					//Malha M4: dx=dy=0.0125
					Nx=1280;
					Ny=480;
				break;
				default:
					printf("Selecionar a opcao de malha corretamente\n");
					exit(1);
				break;
			}
		break;
		case 4://Contracao 5:1
			//8L x 8L
			x0=0.0;
			xm=8.0;
			xf=16.0; 
			y0=0.0;
		        //ym=4.0;	
			ymDown=4.0;
			ymUp=6.0;
			yf=10.0;
			ym = 0.5*(yf+y0);
			aInlet = -3.0/10.0;
			switch(Malha){
				case 0:
					//Malha M1: dx=dy=0.1
					Nx=160; 
					Ny=100;
				break;
				case 1:
					//Malha M2: dx=dy=0.05
					Nx=320;
					Ny=200;
				break;
				case 2:
					//Malha M3: dx=dy=0.025
					Nx=640;
					Ny=400;
				break;
				case 3:
					//Malha M4: dx=dy=0.0125
					Nx=1280;
					Ny=800;
				break;
				default:
					printf("Selecionar a opcao de malha corretamente\n");
					exit(1);
				break;
			}

		break;
		default:
			printf("Selecionar a opcao de geometria correta\n");
			exit(1);
		break;
	}

	double *x, *y, *t; 
	double **uOld, **vOld, **psiOld, **omegaOld;
	double **uNew, **vNew, **psiNew, **omegaNew;
	double **pOld, **pNew;
	double **txxOld, **txyOld, **tyyOld, **txxNew, **txyNew, **tyyNew;//elastic stress (13/03/2026)
	double **dudx, **dudy, **dvdx, **dvdy;
	double *ErroVec, *OmegaVec, tol=1.0e-06, NormErroVec, NormOmegaVec, ErroRel;
	int **VecInd;
	double c, ConvU, ConvV, ViscX, ViscY;
	double dx, dy, dt;
	int i, j, n, k, p, dim;
	int i0, im, ie;		//i0: inicial; im: meio; ie: final
	int j0, jm1, jm2, je;
	indices *vec;
	double Re, beta2, aux, kappa_slip, kappa_Not;
       	double Wi, beta, rho, U=1.0, L=1.0, mu0, mus, mup, lambda, alpha, epsilon, Visc;	//viscoelastic flow parameters (13/03/2026)
	int slip;
	int LinearSystem;

	//Parameters:
	//_______________________________________
	//L scale (characteristic length scale):
	L = 1.0;
	//U scale (velocity length scale):
	U = 1.0;
	//density (ou massa especifica):
	rho = 1.0;
	//Viscosidade do solvente:
	mus = 0.5;
	//Viscosidade do polimero:
	mup = 0.5;
	//Viscosidade dinamica total (mu0):
	mu0 = mus + mup;
	//Razao entre viscosidade do solvente pela viscosidade total (beta):
	beta = mus/mu0;
	//Tempo de relaxacao (lambda):
	lambda = 1.0;
	//Reynolds number:
	Re = (U*L*rho)/mu0;//10.0;
	//Weissenberg number:
	Wi = lambda * (U/L) ;
	//PTT viscoelastic model (setar zero para outros modelos):
	epsilon = 0.25;
	//Giesekus viscoelastic model (setar zero para outros modelos):
	alpha = 0.0;
	//Viscoelastic flow:
	//Ligado	:	1.0
	//Desligado	:	0.0
	Visc = 1.0;
	//--------------------------------------
	//_______________________________
	//Gauss-Seidel:			0
	//Conjugate Gradient Method:	1
	LinearSystem = 1;
	//Espacamento da Malha:
	dx = (xf-x0)/Nx;
	dy = (yf-y0)/Ny;
	dt = (tf-t0)/Nt;
	//_____________________
	//Condicao Navier-slip:
	//_____________________
	//Desligado (No-slip):				0
	//Ligado (Navier-slip):				1
	//Slip modificado: kappa_slip = kappa_Not * r:	2
	slip = 0;
	//____________________
	//Slip coeficient:
	kappa_Not = 100.0; 
	kappa_slip = 1.0;
	//----------------------------
	//Viscoelastic parameters:
	//
	printf("Parametros da Simulacao:\n");
	printf("x0=%f, xm=%f, xf=%f\n", x0, xm, xf);
	printf("y0=%f, ymDown=%f, ymUp=%f, yf=%f\n", y0, ymDown, ymUp, yf);
	printf("t0=%f, tf=%f\n", t0, tf);
	printf("dx=%f, dy=%f, dt=%f\n\n", dx, dy, dt);//getchar();

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
	txxOld = (double **)malloc(sizeof(double *)*(Nx+1));
	txyOld = (double **)malloc(sizeof(double *)*(Nx+1));
	tyyOld = (double **)malloc(sizeof(double *)*(Nx+1));

	dudx = (double **)malloc(sizeof(double *)*(Nx+1));
	dudy = (double **)malloc(sizeof(double *)*(Nx+1));
	dvdx = (double **)malloc(sizeof(double *)*(Nx+1));
	dvdy = (double **)malloc(sizeof(double *)*(Nx+1));

	//--------------------------------------------------
	uNew = (double **)malloc(sizeof(double *)*(Nx+1));
	vNew = (double **)malloc(sizeof(double *)*(Nx+1));
	psiNew = (double **)malloc(sizeof(double *)*(Nx+1));
	omegaNew = (double **)malloc(sizeof(double *)*(Nx+1));
	pNew = (double **)malloc(sizeof(double *)*(Nx+1));
	txxNew = (double **)malloc(sizeof(double *)*(Nx+1));
	txyNew = (double **)malloc(sizeof(double *)*(Nx+1));
	tyyNew = (double **)malloc(sizeof(double *)*(Nx+1));

	VecInd = (int **)malloc(sizeof(int *)*(Nx + 1));
	for (i=0; i<=Nx; i++){
		uOld[i] = (double *)malloc(sizeof(double)*(Ny+1));
		vOld[i] = (double *)malloc(sizeof(double)*(Ny+1));
		psiOld[i] = (double *)malloc(sizeof(double)*(Ny+1));
		omegaOld[i] = (double *)malloc(sizeof(double)*(Ny+1));
		pOld[i] = (double *)malloc(sizeof(double)*(Ny+1));
		txxOld[i] = (double *)malloc(sizeof(double)*(Ny+1));
		txyOld[i] = (double *)malloc(sizeof(double)*(Ny+1));
		tyyOld[i] = (double *)malloc(sizeof(double)*(Ny+1));

		dudx[i] = (double *)malloc(sizeof(double)*(Ny+1));
		dudy[i] = (double *)malloc(sizeof(double)*(Ny+1));
		dvdx[i] = (double *)malloc(sizeof(double)*(Ny+1));
		dvdy[i] = (double *)malloc(sizeof(double)*(Ny+1));

		//--------------------------------------------------
		uNew[i] = (double *)malloc(sizeof(double)*(Ny+1));
		vNew[i] = (double *)malloc(sizeof(double)*(Ny+1));
		psiNew[i] = (double *)malloc(sizeof(double)*(Ny+1));
		omegaNew[i] = (double *)malloc(sizeof(double)*(Ny+1));
		pNew[i] = (double *)malloc(sizeof(double)*(Ny+1));
		txxNew[i] = (double *)malloc(sizeof(double)*(Ny+1));
		txyNew[i] = (double *)malloc(sizeof(double)*(Ny+1));
		tyyNew[i] = (double *)malloc(sizeof(double)*(Ny+1));

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
			txxOld[i][j] = txxNew[i][j] = 0.0;
			txyOld[i][j] = txyNew[i][j] = 0.0;
			tyyOld[i][j] = tyyNew[i][j] = 0.0;
			dudx[i][j] = 0.0;
			dudy[i][j] = 0.0;
			dvdx[i][j] = 0.0;
			dvdy[i][j] = 0.0;
		}
	}
	//printf("\n apos inicializacao Ny=%d\n\n \n", Ny);
	//Boundary conditions:
	//INFLOW:
	i=0;
	//constntant c:
	c = 0.0;
	for (j=0; j<=Ny; j++){

		//printf("j=%d\n", j);
		//uOld[i][j] = uNew[i][j] = -(3.0/32.0)*( y[j] - 4.0 )*y[j];
		//uOld[i][j] = uNew[i][j] = -(3.0/(8.0*(ym*ym)))*( y[j] - yf )*y[j];//Generalizado - 22/04/2024
		//uOld[i][j] = uNew[i][j] = aInlet*(1.0/(ym*ym))*( y[j] - yf )*y[j];//Generalizado - 15/10/2024 (Funcionando)
		uOld[i][j] = uNew[i][j] = ( Umax/((ym-y0)*(ym-yf)) )*(y[j] - y0)*(y[j] - yf);//Generalizado - 10/03/2026 (Funcionando)

		vOld[i][j] = vNew[i][j] = 0.0;
		//psiOld[i][j] = psiNew[i][j] = -(1.0/32.0)*y[j]*y[j]*y[j] + (3.0/16.0)*y[j]*y[j] + c;
		//psiOld[i][j] = psiNew[i][j] = -(3.0/(8.0*(ym*ym)))*( (y[j]/3.0) - (yf/2.0) )*y[j]*y[j] + c;//Generalizado - 22/04/2024 (Funcionando)
		//psiOld[i][j] = psiNew[i][j] = aInlet*(1.0/(ym*ym))*( (y[j]/3.0) - (yf/2.0) )*y[j]*y[j] + c;//Generalizado - 15/10/2024 (Funcionando)
		psiOld[i][j] = psiNew[i][j] = ( Umax/((ym-y0)*(ym-yf)) )*( ((y[j]*y[j]*y[j])/3.0) - (y0 + yf)*((y[j]*y[j])/2.0) + y0*yf*y[j] - ((y0*y0*y0)/3.0) + (y0 + yf)*((y0*y0)/2.0) - y0*yf*y0);//Generalizado - 10/03/2026


		//psiOld[i][j] = psiOld[i+1][j] + dy*uOld[i][j];
		//psiNew[i][j] = psiNew[i+1][j] + dy*uNew[i][j];
		//omegaOld[i][j] = omegaNew[i][j] = (3.0/16.0)*y[j] - (3.0/8.0);
		//omegaOld[i][j] = omegaNew[i][j] = -(3.0/(8.0*(ym*ym)))*( 2.0*y[j] - yf );//Generalizado - 22/04/2024
		//omegaOld[i][j] = omegaNew[i][j] = aInlet*(1.0/(ym*ym))*( 2.0*y[j] - yf );//Generalizado - 15/10/2024 (Funcionando)
		omegaOld[i][j] = omegaNew[i][j] = ( Umax/((ym-y0)*(ym-yf)) )*( 2.0*y[j] - (y0 + yf) );//Generalizado - 10/03/2026 (Funcionando)

		//Oldroyd-B analytical solution:
		//txx = 2*Wi(1-beta)*(dudy)^2
		//txy = (1-beta)*dudy
		//tyy = 0
		dudy[i][j] = (Umax/((ym - y0)*(ym - yf)))*( 2.0*y[j] - (y0 + yf) );
		txxOld[i][j] = txxNew[i][j] = 2.0*Wi*(1.0 - beta)*dudy[i][j]*dudy[i][j];
		txyOld[i][j] = txyNew[i][j] = (1.0 - beta)*dudy[i][j];
		tyyOld[i][j] = tyyNew[i][j] = 0.0;
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
	
	ErroVec = (double *)malloc(sizeof(double)*(dim));
	OmegaVec = (double *)malloc(sizeof(double)*(dim));

	/*	
	for (k=0; k<dim; k++){
		printf("k=%d\n", k);
		printf("i=%d, j=%d\n", vec[k].i, vec[k].j);getchar();
	}
	*/
	

	/*
	//Dados para a resolucao da equacao de Poisson:
	double **A, *s, *b;
	int iter;
	A=(double **)malloc(sizeof(double *)*(dim+1));
	for (k=0; k<(dim+1); k++){
		A[k]=(double *)malloc(sizeof(double)*(dim+1));
	}
	b = (double *)malloc(sizeof(double)*(dim+1));
	s = (double *)malloc(sizeof(double)*(dim+1));
	*/

	//printf("salvar\n");
	//Salvando estado inicial
	//if (n % 100 == 0.0){
	ImprimeArquivoVTK(uNew, vNew, pNew, psiNew, omegaNew, txxNew, txyNew, tyyNew, x0, xf, y0, yf, x, y, Nx, Ny, n);
	//}
	//Parametros da condicao Navier-slip: 11/03/2026
	double alphax, betax, gamma, a0, a1, a2, a3, param1, param2, param3, param4, radius;
	double divTau;
	//printf("depois de salvar\n");
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
		//_____________
		//PAREDE BAIXO:
		j=j0;
		for (i=0; i<=im; i++){

			if (slip == 2){//Navier-slip + modificação sugerida pelo Jonathan 12/03/2026
				radius = sqrt( (x[i] - xm)*(x[i] - xm) + (y[j] - ymDown)*(y[j] - ymDown) );
				if (radius <= 1.0){
					kappa_slip = kappa_Not * radius;
				}else{
					kappa_slip = kappa_Not;
				}
			}
			//Calculo dos pesos para o uso da condicao Navier-slip
			alphax = (kappa_slip/Re)*dy + ((dy*dy)/2.0);
			betax = 2.0*(kappa_slip/Re)*dy + 4.0*((dy*dy)/2.0);
			gamma = 3.0*(kappa_slip/Re)*dy + 9.0*((dy*dy)/2.0);
			//alguns calculos:
			param1 = (2.0/3.0) - (betax/alphax)*(1.0/24.0);
			param2 = (4.0/3.0) - (betax/alphax)*(1.0/6.0);
			param3 = (9.0/2.0) - (gamma/alphax)*(1.0/6.0);
			param4 = (27.0/8.0) - (gamma/alphax)*(1.0/24.0); 
			//Pesos:
			a3 = ( -(1.0/(24.0*alphax)) + (param1/param2)*(1.0/(6.0*alphax)) )/(param4 - (param1/param2)*param3);
			a2 = ( -(1.0/(6.0*alphax)) - param3*a3 )/param2;
			a1 = (1.0/alphax) - (betax/alphax)*a2 - (gamma/alphax)*a3;
			a0 = - a1 - a2 - a3;

			//omegaNew[i][j] = omegaOld[i][j] = (-85.0*psiOld[i][j] + 108.0*psiOld[i][j+1] - 27.0*psiOld[i][j+2] + 4.0*psiOld[i][j+3])/(18.0*dy*dy);
			//omegaOld[i][j] = - (psiOld[i][j] - 2.0*psiOld[i][j+1] + psiOld[i][j+2])/(dy*dy);
			omegaNew[i][j] = omegaOld[i][j] = -(1.0/( 18.0*(dy*dy)) )*( -85.0*psiOld[i][j] + 108.0*psiOld[i][j+1] - 27.0*psiOld[i][j+2] + 4.0*psiOld[i][j+3] );//Goiabinha
			//omegaNew[i][j] = omegaOld[i][j] = - (1.0/( 18.0*(dy*dy)) )*( 2.0*psiOld[i][j] - 5.0*psiOld[i][j+1] + 4.0*psiOld[i][j+2] - 1.0*psiOld[i][j+3] );//O(h^2) 10/03/2026
			if (slip != 0){//Condição Navier-slip para omega: 10/03/2026
				//omegaNew[i][j] = omegaOld[i][j] = -(Re/kappa_slip)*( ( psiOld[i][j+1] - psiOld[i][j])/dy);//Teste 1 - fail
				//omegaNew[i][j] = omegaOld[i][j] = - 2.0*( ( psiOld[i][j+1] - psiOld[i][j])/( ((2.0*dy*kappa_slip)/Re) + dy*dy) );//Teste 2-fail
				omegaNew[i][j] = omegaOld[i][j] = - (a0*psiOld[i][j] + a1*psiOld[i][j+1] + a2*psiOld[i][j+2] + a3*psiOld[i][j+3]);//Teste 3
			}
			//Extrapolacao linear:
			txxOld[i][j] = txxNew[i][j] = 2.0*txxOld[i][j+1] - txxOld[i][j+2];
			txyOld[i][j] = txyNew[i][j] = 2.0*txyOld[i][j+1] - txyOld[i][j+2];
			tyyOld[i][j] = tyyNew[i][j] = 2.0*tyyOld[i][j+1] - tyyOld[i][j+2];
		}
		j=jm1;
		for (i=im; i<=ie; i++){

			if (slip == 2){//Navier-slip + modificação sugerida pelo Jonathan 12/03/2026
				radius = sqrt( (x[i] - xm)*(x[i] - xm) + (y[j] - ymDown)*(y[j] - ymDown) );
				if (radius <= 1.0){
					kappa_slip = kappa_Not * radius;
				}else{
					kappa_slip = kappa_Not;
				}
			}
			//Calculo dos pesos para o uso da condicao Navier-slip
			alphax = (kappa_slip/Re)*dy + ((dy*dy)/2.0);
			betax = 2.0*(kappa_slip/Re)*dy + 4.0*((dy*dy)/2.0);
			gamma = 3.0*(kappa_slip/Re)*dy + 9.0*((dy*dy)/2.0);
			//alguns calculos:
			param1 = (2.0/3.0) - (betax/alphax)*(1.0/24.0);
			param2 = (4.0/3.0) - (betax/alphax)*(1.0/6.0);
			param3 = (9.0/2.0) - (gamma/alphax)*(1.0/6.0);
			param4 = (27.0/8.0) - (gamma/alphax)*(1.0/24.0); 
			//Pesos:
			a3 = ( -(1.0/(24.0*alphax)) + (param1/param2)*(1.0/(6.0*alphax)) )/(param4 - (param1/param2)*param3);
			a2 = ( -(1.0/(6.0*alphax)) - param3*a3 )/param2;
			a1 = (1.0/alphax) - (betax/alphax)*a2 - (gamma/alphax)*a3;
			a0 = - a1 - a2 - a3;

			omegaNew[i][j] = omegaOld[i][j] = -(-85.0*psiOld[i][j] + 108.0*psiOld[i][j+1] - 27.0*psiOld[i][j+2] + 4.0*psiOld[i][j+3])/(18.0*dy*dy);
			//omegaNew[i][j] = omegaOld[i][j] = - ( 2.0*psiOld[i][j] - 5.0*psiOld[i][j+1] + 4.0*psiOld[i][j+2] - 1.0*psiOld[i][j+3])/(dy*dy);//O(h^2) 10/03/2026
			//omegaOld[i][j] = - (psiOld[i][j] - 2.0*psiOld[i][j+1] + psiOld[i][j+2])/(dy*dy);
			if (slip != 0){//Condição Navier-slip para omega: 10/03/2026
				//omegaNew[i][j] = omegaOld[i][j] = -(Re/kappa_slip)*( ( psiOld[i][j+1] - psiOld[i][j])/dy);//Teste 1 - fail
				//omegaNew[i][j] = omegaOld[i][j] = - 2.0*( ( psiOld[i][j+1] - psiOld[i][j])/( ((2.0*dy*kappa_slip)/Re) + dy*dy) );//Teste 2-fail
				omegaNew[i][j] = omegaOld[i][j] = - (a0*psiOld[i][j] + a1*psiOld[i][j+1] + a2*psiOld[i][j+2] + a3*psiOld[i][j+3]);//Teste 3
			}
			//Extrapolacao linear:
			txxOld[i][j] = txxNew[i][j] = 2.0*txxOld[i][j+1] - txxOld[i][j+2];
			txyOld[i][j] = txyNew[i][j] = 2.0*txyOld[i][j+1] - txyOld[i][j+2];
			tyyOld[i][j] = tyyNew[i][j] = 2.0*tyyOld[i][j+1] - tyyOld[i][j+2];
		}
		//printf("PAREDE BAIXO\n");
		//________________
		//PAREDE SUPERIOR:
		j=je;
		for (i=0; i<=im; i++){

			if (slip == 2){//Navier-slip + modificação sugerida pelo Jonathan 12/03/2026
				radius = sqrt( (x[i] - xm)*(x[i] - xm) + (y[j] - ymDown)*(y[j] - ymDown) );
				if (radius <= 1.0){
					kappa_slip = kappa_Not * radius;
				}else{
					kappa_slip = kappa_Not;
				}
			}
			//Calculo dos pesos para o uso da condicao Navier-slip
			alphax = (kappa_slip/Re)*dy + ((dy*dy)/2.0);
			betax = 2.0*(kappa_slip/Re)*dy + 4.0*((dy*dy)/2.0);
			gamma = 3.0*(kappa_slip/Re)*dy + 9.0*((dy*dy)/2.0);
			//alguns calculos:
			param1 = (2.0/3.0) - (betax/alphax)*(1.0/24.0);
			param2 = (4.0/3.0) - (betax/alphax)*(1.0/6.0);
			param3 = (9.0/2.0) - (gamma/alphax)*(1.0/6.0);
			param4 = (27.0/8.0) - (gamma/alphax)*(1.0/24.0); 
			//Pesos:
			a3 = ( -(1.0/(24.0*alphax)) + (param1/param2)*(1.0/(6.0*alphax)) )/(param4 - (param1/param2)*param3);
			a2 = ( -(1.0/(6.0*alphax)) - param3*a3 )/param2;
			a1 = (1.0/alphax) - (betax/alphax)*a2 - (gamma/alphax)*a3;
			a0 = - a1 - a2 - a3;

			omegaNew[i][j] = omegaOld[i][j] = -(-85.0*psiOld[i][j] + 108.0*psiOld[i][j-1] - 27.0*psiOld[i][j-2] + 4.0*psiOld[i][j-3])/(18.0*dy*dy);
			//omegaOld[i][j] = - (psiOld[i][j] - 2.0*psiOld[i][j-1] + psiOld[i][j-2])/(dy*dy);
			//omegaNew[i][j] = omegaOld[i][j] = - ( 2.0*psiOld[i][j] - 5.0*psiOld[i][j-1] + 4.0*psiOld[i][j-2] - 1.0*psiOld[i][j-3])/(dy*dy);//O(h^2) 10/03/2026
			if (slip != 0){//Condição Navier-slip para omega: 10/03/2026
				//omegaNew[i][j] = omegaOld[i][j] = (Re/kappa_slip)*( ( psiOld[i][j] - psiOld[i][j-1])/dy);//Teste 1 - fail
				//omegaNew[i][j] = omegaOld[i][j] = - 2.0*( ( psiOld[i][j] - psiOld[i][j-1])/( ((2.0*dy*kappa_slip)/Re) - dy*dy) );//Teste 2 - fail
				omegaNew[i][j] = omegaOld[i][j] = - (a0*psiOld[i][j] + a1*psiOld[i][j-1] + a2*psiOld[i][j-2] + a3*psiOld[i][j-3]);//Teste 3
			}
			//Extrapolacao linear:
			txxOld[i][j] = txxNew[i][j] = 2.0*txxOld[i][j-1] - txxOld[i][j-2];
			txyOld[i][j] = txyNew[i][j] = 2.0*txyOld[i][j-1] - txyOld[i][j-2];
			tyyOld[i][j] = tyyNew[i][j] = 2.0*tyyOld[i][j-1] - tyyOld[i][j-2];
		}
		j=jm2;
		for (i=im; i<=ie; i++){

			if (slip == 2){//Navier-slip + modificação sugerida pelo Jonathan 12/03/2026
				radius = sqrt( (x[i] - xm)*(x[i] - xm) + (y[j] - ymDown)*(y[j] - ymDown) );
				if (radius <= 1.0){
					kappa_slip = kappa_Not * radius;
				}else{
					kappa_slip = kappa_Not;
				}
			}
			//Calculo dos pesos para o uso da condicao Navier-slip
			alphax = (kappa_slip/Re)*dy + ((dy*dy)/2.0);
			betax = 2.0*(kappa_slip/Re)*dy + 4.0*((dy*dy)/2.0);
			gamma = 3.0*(kappa_slip/Re)*dy + 9.0*((dy*dy)/2.0);
			//alguns calculos:
			param1 = (2.0/3.0) - (betax/alphax)*(1.0/24.0);
			param2 = (4.0/3.0) - (betax/alphax)*(1.0/6.0);
			param3 = (9.0/2.0) - (gamma/alphax)*(1.0/6.0);
			param4 = (27.0/8.0) - (gamma/alphax)*(1.0/24.0); 
			//Pesos:
			a3 = ( -(1.0/(24.0*alphax)) + (param1/param2)*(1.0/(6.0*alphax)) )/(param4 - (param1/param2)*param3);
			a2 = ( -(1.0/(6.0*alphax)) - param3*a3 )/param2;
			a1 = (1.0/alphax) - (betax/alphax)*a2 - (gamma/alphax)*a3;
			a0 = - a1 - a2 - a3;

			omegaNew[i][j] = omegaOld[i][j] = -(-85.0*psiOld[i][j] + 108.0*psiOld[i][j-1] - 27.0*psiOld[i][j-2] + 4.0*psiOld[i][j-3])/(18.0*dy*dy);
			//omegaOld[i][j] = - (psiOld[i][j] - 2.0*psiOld[i][j-1] + psiOld[i][j-2])/(dy*dy);
			//omegaNew[i][j] = omegaOld[i][j] = - ( 2.0*psiOld[i][j] - 5.0*psiOld[i][j-1] + 4.0*psiOld[i][j-2] - 1.0*psiOld[i][j-3])/(dy*dy);//O(h^2) 10/03/2026
			if (slip != 0){//Condição Navier-slip para omega: 10/03/2026
				//omegaNew[i][j] = omegaOld[i][j] = (Re/kappa_slip)*( ( psiOld[i][j] - psiOld[i][j-1])/dy);//Teste 1 - fail
				//omegaNew[i][j] = omegaOld[i][j] = - 2.0*( ( psiOld[i][j] - psiOld[i][j-1])/( ((2.0*dy*kappa_slip)/Re) - dy*dy) );//Teste 2
				omegaNew[i][j] = omegaOld[i][j] = - (a0*psiOld[i][j] + a1*psiOld[i][j-1] + a2*psiOld[i][j-2] + a3*psiOld[i][j-3]);//Teste 3
			}
			//Extrapolacao linear:
			txxOld[i][j] = txxNew[i][j] = 2.0*txxOld[i][j-1] - txxOld[i][j-2];
			txyOld[i][j] = txyNew[i][j] = 2.0*txyOld[i][j-1] - txyOld[i][j-2];
			tyyOld[i][j] = tyyNew[i][j] = 2.0*tyyOld[i][j-1] - tyyOld[i][j-2];
		}
		//printf("PAREDE SUPERIOR\n");
		//_______________
		//PAREDE DIREITA:
		i=im;
		//Parte de baixo:
		for (j=0; j<=jm1; j++){

			if (slip == 2){//Navier-slip + modificação sugerida pelo Jonathan 12/03/2026
				radius = sqrt( (x[i] - xm)*(x[i] - xm) + (y[j] - ymDown)*(y[j] - ymDown) );
				if (radius <= 1.0){
					kappa_slip = kappa_Not * radius;
				}else{
					kappa_slip = kappa_Not;
				}
			}
			//Calculo dos pesos para o uso da condicao Navier-slip
			alphax = (kappa_slip/Re)*dx + ((dx*dx)/2.0);
			//printf("(kappa_slip/Re)*dx=%f, ((dx*dx)/2.0)=%f\n", (kappa_slip/Re)*dx, ((dx*dx)/2.0));getchar();
			betax = 2.0*(kappa_slip/Re)*dx + 4.0*((dx*dx)/2.0);
			gamma = 3.0*(kappa_slip/Re)*dx + 9.0*((dx*dx)/2.0);
			//printf("alphax=%f, betax=%f, gamma=%f\n", alphax, betax, gamma);getchar();
			//alguns calculos:
			param1 = (2.0/3.0) - (betax/alphax)*(1.0/24.0);
			param2 = (4.0/3.0) - (betax/alphax)*(1.0/6.0);
			param3 = (9.0/2.0) - (gamma/alphax)*(1.0/6.0);
			param4 = (27.0/8.0) - (gamma/alphax)*(1.0/24.0); 
			//Pesos:
			a3 = ( -(1.0/(24.0*alphax)) + (param1/param2)*(1.0/(6.0*alphax)) )/(param4 - (param1/param2)*param3);
			a2 = ( -(1.0/(6.0*alphax)) - param3*a3 )/param2;
			a1 = (1.0/alphax) - (betax/alphax)*a2 - (gamma/alphax)*a3;
			a0 = - a1 - a2 - a3;

			//printf("im=%d, j=%d\n", i, j);
			omegaNew[i][j] = omegaOld[i][j] = -(-85.0*psiOld[i][j] + 108.0*psiOld[i-1][j] - 27.0*psiOld[i-2][j] + 4.0*psiOld[i-3][j])/(18.0*dx*dx);
			//omegaOld[i][j] = -(psiOld[i][j] - 2.0*psiOld[i-1][j] + psiOld[i-2][j])/(dx*dx);
			//omegaNew[i][j] = omegaOld[i][j] = - ( 2.0*psiOld[i][j] - 5.0*psiOld[i-1][j] + 4.0*psiOld[i-2][j] - 1.0*psiOld[i-3][j] )/(dx*dx);//O(h^2) 10/03/2026
			//printf("omegaOld[%d][%d]=%f\n", i, j, omegaOld[i][j]);
			if (slip != 0){//Condição Navier-slip para omega: 10/03/2026
				//omegaNew[i][j] = omegaOld[i][j] = (Re/kappa_slip)*( ( psiOld[i][j] - psiOld[i-1][j])/dx);//Teste 1 - fail
				//omegaNew[i][j] = omegaOld[i][j] = -2.0*( ( psiOld[i-1][j] - psiOld[i][j])/( ((2.0*dx*kappa_slip)/Re) + dx*dx) );//Teste 2-fail
				omegaNew[i][j] = omegaOld[i][j] = - (a0*psiOld[i][j] + a1*psiOld[i-1][j] + a2*psiOld[i-2][j] + a3*psiOld[i-3][j]);//Teste 3
			}
			//Extrapolacao linear:
			txxOld[i][j] = txxNew[i][j] = 2.0*txxOld[i-1][j] - txxOld[i-2][j];
			txyOld[i][j] = txyNew[i][j] = 2.0*txyOld[i-1][j] - txyOld[i-2][j];
			tyyOld[i][j] = tyyNew[i][j] = 2.0*tyyOld[i-1][j] - tyyOld[i-2][j];
		}
		for (j=jm2; j<=je; j++){

			if (slip == 2){//Navier-slip + modificação sugerida pelo Jonathan 12/03/2026
				radius = sqrt( (x[i] - xm)*(x[i] - xm) + (y[j] - ymDown)*(y[j] - ymDown) );
				if (radius <= 1.0){
					kappa_slip = kappa_Not * radius;
				}else{
					kappa_slip = kappa_Not;
				}
			}
			//Calculo dos pesos para o uso da condicao Navier-slip
			alphax = (kappa_slip/Re)*dx + ((dx*dx)/2.0);
			//printf("(kappa_slip/Re)*dx=%f, ((dx*dx)/2.0)=%f\n", (kappa_slip/Re)*dx, ((dx*dx)/2.0));getchar();
			betax = 2.0*(kappa_slip/Re)*dx + 4.0*((dx*dx)/2.0);
			gamma = 3.0*(kappa_slip/Re)*dx + 9.0*((dx*dx)/2.0);
			//printf("alphax=%f, betax=%f, gamma=%f\n", alphax, betax, gamma);getchar();
			//alguns calculos:
			param1 = (2.0/3.0) - (betax/alphax)*(1.0/24.0);
			param2 = (4.0/3.0) - (betax/alphax)*(1.0/6.0);
			param3 = (9.0/2.0) - (gamma/alphax)*(1.0/6.0);
			param4 = (27.0/8.0) - (gamma/alphax)*(1.0/24.0); 
			//Pesos:
			a3 = ( -(1.0/(24.0*alphax)) + (param1/param2)*(1.0/(6.0*alphax)) )/(param4 - (param1/param2)*param3);
			a2 = ( -(1.0/(6.0*alphax)) - param3*a3 )/param2;
			a1 = (1.0/alphax) - (betax/alphax)*a2 - (gamma/alphax)*a3;
			a0 = - a1 - a2 - a3;

			omegaNew[i][j] = omegaOld[i][j] = -(-85.0*psiOld[i][j] + 108.0*psiOld[i-1][j] - 27.0*psiOld[i-2][j] + 4.0*psiOld[i-3][j])/(18.0*dx*dx);
			//omegaNew[i][j] = omegaOld[i][j] = - ( 2.0*psiOld[i][j] - 5.0*psiOld[i-1][j] + 4.0*psiOld[i-2][j] - 1.0*psiOld[i-3][j] )/(dx*dx);//O(h^2) 10/03/2026
			//omegaOld[i][j] = -(psiOld[i][j] - 2.0*psiOld[i-1][j] + psiOld[i-2][j])/(dx*dx);
			if (slip != 0){//Condição Navier-slip para omega: 10/03/2026
				//omegaNew[i][j] = omegaOld[i][j] = (Re/kappa_slip)*( ( psiOld[i][j] - psiOld[i-1][j])/dx);//Teste 1 - fail
				//omegaNew[i][j] = omegaOld[i][j] = -2.0*( ( psiOld[i-1][j] - psiOld[i][j])/( ((2.0*dx*kappa_slip)/Re) + dx*dx) );//Teste 2-fail
				omegaNew[i][j] = omegaOld[i][j] = - (a0*psiOld[i][j] + a1*psiOld[i-1][j] + a2*psiOld[i-2][j] + a3*psiOld[i-3][j]);//Teste 3
			}
			//Extrapolacao linear:
			txxOld[i][j] = txxNew[i][j] = 2.0*txxOld[i-1][j] - txxOld[i-2][j];
			txyOld[i][j] = txyNew[i][j] = 2.0*txyOld[i-1][j] - txyOld[i-2][j];
			tyyOld[i][j] = tyyNew[i][j] = 2.0*tyyOld[i-1][j] - tyyOld[i-2][j];
		}
		//printf("PAREDE DIREITA\n");
		//________
		//OUTFLOW:
		i=ie;
		for (j=jm1; j<=jm2; j++){
			omegaNew[i][j] = omegaNew[i-1][j];
			psiNew[i][j] = psiNew[i-1][j];

			omegaOld[i][j] = omegaOld[i-1][j];
			psiOld[i][j] = psiOld[i-1][j];

			txxNew[i][j] = txxNew[i-1][j];
			txyNew[i][j] = txyNew[i-1][j];
			tyyNew[i][j] = tyyNew[i-1][j];
			txxOld[i][j] = txxOld[i-1][j];
			txyOld[i][j] = txyOld[i-1][j];
			tyyOld[i][j] = tyyOld[i-1][j];
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
					ConvU = -uOld[i][j]*(omegaOld[i+1][j] - omegaOld[i-1][j])/(2.0*dx);//DiferenÃ§a centrada			 
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

				divTau = ((txyOld[i+1][j] -2.0*txyOld[i][j] + txyOld[i-1][j])/(dx*dx)) + ((tyyOld[i+1][j+1] - tyyOld[i-1][j+1] - tyyOld[i+1][j-1] + tyyOld[i-1][j-1])/(4.0*dx*dy)) - ((txxOld[i+1][j+1] - txxOld[i-1][j+1] - txxOld[i+1][j-1] + txxOld[i-1][j-1])/(4.0*dx*dy)) - ((txyOld[i][j+1] - 2.0*txyOld[i][j] + txyOld[i][j-1])/(dy*dy));

				//Euler Explicito - omegaNew
				//omegaNew[i][j] = omegaOld[i][j] + dt*( ConvU + ConvV + (1.0/Re)*( ViscX + ViscY ) );
				omegaNew[i][j] = omegaOld[i][j] + dt*( -uOld[i][j]*(omegaOld[i+1][j] - omegaOld[i-1][j])/(2.0*dx) - vOld[i][j]*(omegaOld[i][j+1] - omegaOld[i][j-1])/(2.0*dy) + (beta/Re)*( (omegaOld[i+1][j] - 2.0*omegaOld[i][j] + omegaOld[i-1][j])/(dx*dx) + (omegaOld[i][j+1] - 2.0*omegaOld[i][j] + omegaOld[i][j-1])/(dy*dy) ) + (1.0/Re)*divTau*Visc );
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

				divTau = ((txyOld[i+1][j] -2.0*txyOld[i][j] + txyOld[i-1][j])/(dx*dx)) + ((tyyOld[i+1][j+1] - tyyOld[i-1][j+1] - tyyOld[i+1][j-1] + tyyOld[i-1][j-1])/(4.0*dx*dy)) - ((txxOld[i+1][j+1] - txxOld[i-1][j+1] - txxOld[i+1][j-1] + txxOld[i-1][j-1])/(4.0*dx*dy)) - ((txyOld[i][j+1] - 2.0*txyOld[i][j] + txyOld[i][j-1])/(dy*dy));

				//Euler Explicito - omegaNew
				//omegaNew[i][j] = omegaOld[i][j] + dt*( ConvU + ConvV + (1.0/Re)*( ViscX + ViscY ) );
				omegaNew[i][j] = omegaOld[i][j] + dt*( -uOld[i][j]*(omegaOld[i+1][j] - omegaOld[i-1][j])/(2.0*dx) - vOld[i][j]*(omegaOld[i][j+1] - omegaOld[i][j-1])/(2.0*dy) + (beta/Re)*( (omegaOld[i+1][j] - 2.0*omegaOld[i][j] + omegaOld[i-1][j])/(dx*dx) + (omegaOld[i][j+1] - 2.0*omegaOld[i][j] + omegaOld[i][j-1])/(dy*dy) ) + (1.0/Re)*divTau*Visc );
			}
		}

		//Calculo do erro relativo:
		for (p=0; p<dim; p++){
			//indices:
			i = vec[p].i;
			j = vec[p].j;

			ErroVec[p] = omegaNew[i][j] - omegaOld[i][j];
			OmegaVec[p] = omegaNew[i][j];
		}

		//NormErroVec = sqrt(InnerProduct(ErroVec, ErroVec, dim));
		//NormOmegaVec = sqrt(InnerProduct(OmegaVec, OmegaVec, dim));

		NormErroVec = Norm(ErroVec, dim);
		NormOmegaVec = Norm(OmegaVec, dim);

		if (NormOmegaVec != 0.0){

			ErroRel = NormErroVec/NormOmegaVec;
		}else{
			ErroRel = 1.0;
		}

		//################################################################
		//Sol. Eq. Poisson:
		//################################################################
		//_______________________________________
		//Atualiza condicao de contorno para psi:		
		//PAREDE VERTICAL ESQUERDA (INFLOW):
		i=i0;
		for (j=0; j<=je; j++){
			//psiOld[i][j] = 2.0*psiOld[i+1][j]  - psiOld[i+2][j];
		}
		
		//PAREDE HORIZONTAL SUPERIOR:
		j=je;
		for (i=1; i<=im; i++){
			//psiOld[i][j] = psiNew[i][j] = -(3.0/(8.0*(ym*ym)))*( (yf/3.0) - (yf/2.0) )*yf*yf + c;
 			//psiOld[i][j] = psiNew[i][j] = aInlet*(1.0/(ym*ym))*( (yf/3.0) - (yf/2.0) )*yf*yf + c;//Generalizado - 15/10/2024
			//psiOld[i][j] = 2.0*psiOld[i][j-1] - 1.0*psiOld[i][j-2];
			//psiNew[i][j] = 2.0*psiNew[i][j-1] - 1.0*psiNew[i][j-2];
			psiOld[i][j] = psiNew[i][j] = (Umax/((ym-y0)*(ym-yf)))*((-yf*yf*yf + y0*y0*y0 + 3.0*yf*yf*y0 - 3.0*yf*y0*y0)/6.0);//Generalizado - 10/03/2026
			//if (slip !=0){
				//psiOld[i][j] = psiNew[i][j] = psiOld[i][j-1] + dy*(kappa_slip/Re)*omegaNew[i][j]; 
			//}

		}
		j=jm2;
		for (i=im; i<=ie; i++){
			//psiOld[i][j] = psiNew[i][j] = -(3.0/(8.0*(ym*ym)))*( (yf/3.0) - (yf/2.0) )*yf*yf + c;
			//psiOld[i][j] = psiNew[i][j] = aInlet*(1.0/(ym*ym))*( (yf/3.0) - (yf/2.0) )*yf*yf + c;
			//psiOld[i][j] = psiOld[i][j-1];
			//psiNew[i][j] = psiNew[i][j-1];
			psiOld[i][j] = psiNew[i][j] = (Umax/((ym-y0)*(ym-yf)))*((-yf*yf*yf + y0*y0*y0 + 3.0*yf*yf*y0 - 3.0*yf*y0*y0)/6.0);//Generalizado - 10/03/2026
			//if (slip !=0){
				//psiOld[i][j] = psiNew[i][j] = psiOld[i][j-1] + dy*(kappa_slip/Re)*omegaNew[i][j]; 
			//}

		}
		//PAREDE HORIZONTAL INFERIOR:
		j=j0;
		for (i=1; i<=im; i++){
			//psiOld[i][j] = psiNew[i][j] = -(3.0/(8.0*(ym*ym)))*( (y0/3.0) - (y0/2.0) )*y0*y0 + c;
			//psiOld[i][j] = psiNew[i][j] = aInlet*(1.0/(ym*ym))*( (y0/3.0) - (y0/2.0) )*y0*y0 + c;//
			//psiOld[i][j] = psiOld[i][j+1];
			//psiNew[i][j] = psiNew[i][j+1];
			psiOld[i][j] = psiNew[i][j] = 0.0;//Generalizado - 10/03/2026
			//if (slip !=0){
				//psiOld[i][j] = psiNew[i][j] = psiOld[i][j+1] + dy*(kappa_slip/Re)*omegaNew[i][j]; 
			//}

		}
		j=jm1;
		for (i=im; i<=ie; i++){
			//psiOld[i][j] = psiNew[i][j] = -(3.0/(8.0*(ym*ym)))*( (y0/3.0) - (y0/2.0) )*y0*y0 + c;
			//psiOld[i][j] = psiNew[i][j] = aInlet*(1.0/(ym*ym))*( (y0/3.0) - (y0/2.0) )*y0*y0 + c;//
			//psiOld[i][j] = psiOld[i][j+1];
			//psiNew[i][j] = psiNew[i][j+1];
			psiOld[i][j] = psiNew[i][j] = 0.0;//Generalizado - 10/03/2026
			//if (slip !=0){
				//psiOld[i][j] = psiNew[i][j] = psiOld[i][j+1] + dy*(kappa_slip/Re)*omegaNew[i][j]; 
			//}
		}
		//PAREDE VERTICAL DIREITA:
		i=im;
		for (j=0; j<=jm1; j++){//baixo
			//psiOld[i][j] = psiNew[i][j] = -(3.0/(8.0*(ym*ym)))*( (y0/3.0) - (y0/2.0) )*y0*y0 + c;
			//psiOld[i][j] = psiNew[i][j] = aInlet*(1.0/(ym*ym))*( (y0/3.0) - (y0/2.0) )*y0*y0 + c;//
			//psiOld[i][j] = psiOld[i-1][j];
			//psiNew[i][j] = psiNew[i-1][j];
			psiOld[i][j] = psiNew[i][j] = 0.0;//Generalizado - 10/03/2026
			//if (slip !=0){
				//psiOld[i][j] = psiNew[i][j] = psiOld[i-1][j] + dx*(kappa_slip/Re)*omegaNew[i][j]; 
			//}

		}
		for (j=jm2; j<=je; j++){//cima
			//psiOld[i][j] = psiNew[i][j] = -(3.0/(8.0*(ym*ym)))*( (yf/3.0) - (yf/2.0) )*yf*yf + c;
			//psiOld[i][j] = psiNew[i][j] = aInlet*(1.0/(ym*ym))*( (yf/3.0) - (yf/2.0) )*yf*yf + c;//
			//psiOld[i][j] = psiOld[i-1][j];
			//psiNew[i][j] = psiNew[i-1][j];
			psiOld[i][j] = psiNew[i][j] = (Umax/((ym-y0)*(ym-yf)))*((-yf*yf*yf + y0*y0*y0 + 3.0*yf*yf*y0 - 3.0*yf*y0*y0)/6.0);//Generalizado - 10/03/2026
			//if (slip !=0){
				//psiOld[i][j] = psiNew[i][j] = psiOld[i-1][j] + dx*(kappa_slip/Re)*omegaNew[i][j]; 
			//}
		}
		//PAREDE VERTICAL DIREITA - OUTFLOW: Neumann Homogeneo
		i=ie;
		for (j=jm1; j<=jm2; j++){
			psiOld[i][j] = psiOld[i-1][j];
			psiNew[i][j] = psiNew[i-1][j];
		}
		//________________________________________________________________
	
		/*
		//Inicializa matriz dos coeficientes A, vetor b e chute inicial (tbm soluÃ§Ã£o) s:
		for (i=0; i<dim; i++){
			b[i] = 0.0;
			s[i] = 0.0;
			printf("s[%d] = %f\n", i, s[i]);
			for (j=0; j<dim; j++){
				A[i][j] = 0.0;
				//printf("A[%d][%d]=%f\n", i, j, A[i][j]);
			}
		}
		*/

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
		//Solucao do sistema linear:
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
				printf("Escolha um metodo vÃ¡lido para a resoluÃ§Ã£o do sistema linear\n");
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

			//psiNew[i][j] = aux * ( psiOld[i - 1][j] + psiOld[i + 1][j] + (beta2) * ( psiOld[i][j + 1] + psiOld[i][j - 1] ) - (dx*dx) * omegaNew[i][j]);
			psiNew[i][j] = aux * ( psiOld[i - 1][j] + psiOld[i + 1][j] + (beta2) * ( psiOld[i][j + 1] + psiOld[i][j - 1] ) + (dx*dx) * omegaNew[i][j]);
		}
		
		//_______________________________________
		//Copia vetor s para a matriz psiNew:
		for (p=0; p<dim; p++){
			//printf("%f\n", s[p]);getchar();
			i=vec[p].i;
			j=vec[p].j;
			//psiNew[i][j] = s[p];//ESSA ATUALIZACAO DEVE SER DESCOMENTADA CASO QUEIRA RESOLVER POR SISTEMA LINEAR

			//Calculo das velocidades:
			//uNew[i][j] = (psiNew[i][j+1] - psiNew[i][j-1])/(2.0*dy);
			//vNew[i][j] = -(psiNew[i+1][j] - psiNew[i-1][j])/(2.0*dx);
			
		  	uNew[i][j] =   (psiNew[i][j+1] - psiNew[i][j-1])/(2.0*dy);
                        vNew[i][j] = - (psiNew[i+1][j] - psiNew[i-1][j])/(2.0*dx);

			//Copia soluÃ§Ã£o nova para velha:
			omegaOld[i][j] = omegaNew[i][j];
			psiOld[i][j] = psiNew[i][j];
			uOld[i][j] = uNew[i][j];
			vOld[i][j] = vNew[i][j];	
			pOld[i][j] = pNew[i][j];

		}
		/*
		for (i=1; i<ie; i++){
			for (j=1; j<je; j++){
				//Copia soluÃ§Ã£o nova para velha:
				omegaOld[i][j] = omegaNew[i][j];
				psiOld[i][j] = psiNew[i][j];
				uOld[i][j] = uNew[i][j];
				vOld[i][j] = vNew[i][j];	
				pOld[i][j] = pNew[i][j];
			}
		}
		*/

		//Equacao constitutiva:
		double convUtxx, convVtxx, deformtxx, PTTxx, Giesekusxx, trTau;
		double convUtxy, convVtxy, deformtxy, PTTxy, Giesekusxy;
		double convUtyy, convVtyy, deformtyy, PTTyy, Giesekusyy;
		trTau = txxOld[i][j] + tyyOld[i][j];

		for (p=0; p<dim; p++){
			//printf("%f\n", s[p]);getchar();
			i=vec[p].i;
			j=vec[p].j;

			dudx[i][j] = (uNew[i+1][j]-uNew[i-1][j])/(2.0*dx);
			dudy[i][j] = (uNew[i][j+1]-uNew[i][j-1])/(2.0*dy);
			dvdx[i][j] = (vNew[i+1][j]-vNew[i-1][j])/(2.0*dx);
			dvdy[i][j] = (vNew[i][j+1]-vNew[i][j-1])/(2.0*dy);

			convUtxx = uNew[i][j]*( (txxOld[i+1][j]-txxOld[i-1][j])/(2.0*dx) );
			convVtxx = vNew[i][j]*( (txxOld[i][j+1]-txxOld[i][j-1])/(2.0*dy) );
			deformtxx = 2.0*(dudx[i][j]*txxOld[i][j] + dudy[i][j]*txyOld[i][j]);
			PTTxx =  (epsilon/(1.0-beta))*trTau*txxOld[i][j];
			Giesekusxx = (alpha/(1.0-beta))*(txxOld[i][j]*txxOld[i][j] + txyOld[i][j]*txyOld[i][j]);
			txxNew[i][j] = txxOld[i][j] + dt*( - convUtxx - convVtxx + deformtxx - (1.0/Wi)*txxOld[i][j] + 2.0*((1.0-beta)/Wi)*dudx[i][j] - PTTxx - Giesekusxx);

			convUtxy = uNew[i][j]*( (txyOld[i+1][j]-txyOld[i-1][j])/(2.0*dx) );
			convVtxy = vNew[i][j]*( (txyOld[i][j+1]-txyOld[i][j-1])/(2.0*dy) );
			deformtxy = dvdx[i][j]*txxOld[i][j] + dudy[i][j]*tyyOld[i][j];
			PTTxy =  (epsilon/(1.0-beta))  *trTau*txyOld[i][j];
			Giesekusxy = (alpha/(1.0-beta))*trTau*txyOld[i][j];
			txyNew[i][j] = txyOld[i][j] + dt*( - convUtxy - convVtxy + deformtxy - (1.0/Wi)*txyOld[i][j] + ((1.0-beta)/Wi)*(dudy[i][j] + dvdx[i][j]) - PTTxy - Giesekusxy);

			convUtyy = uNew[i][j]*( (tyyOld[i+1][j]-tyyOld[i-1][j])/(2.0*dx) );
			convVtyy = vNew[i][j]*( (tyyOld[i][j+1]-tyyOld[i][j-1])/(2.0*dy) );
			deformtyy = 2.0*(dvdx[i][j]*txyOld[i][j] + dvdy[i][j]*tyyOld[i][j]);
			PTTyy =  (epsilon/(1.0-beta))*trTau*tyyOld[i][j];
			Giesekusyy = (alpha/(1.0-beta))*(tyyOld[i][j]*tyyOld[i][j] + txyOld[i][j]*txyOld[i][j]);
			tyyNew[i][j] = tyyOld[i][j] + dt*( - convUtyy - convVtyy + deformtyy - (1.0/Wi)*tyyOld[i][j] + 2.0*((1.0-beta)/Wi)*dvdy[i][j] - PTTyy - Giesekusyy);
		}

		//Atualiza tensoes:
		for (p=0; p<dim; p++){
			//printf("%f\n", s[p]);getchar();
			i=vec[p].i;
			j=vec[p].j;

			txxOld[i][j] = txxNew[i][j];
			txyOld[i][j] = txyNew[i][j];
			tyyOld[i][j] = tyyNew[i][j];
		}

		printf("n=%d, ErroRel=%e\n", n, ErroRel);
		//Salvando valores:
		if (n % (Nt/10) == 0.0){
			ImprimeArquivoVTK(uNew, vNew, pNew, psiNew, omegaNew, txxNew, txyNew, tyyNew, x0, xf, y0, yf, x, y, Nx, Ny, n);
		}

	}//end while

	/*
	//Vetores da equaÃ§Ã£o de Poisson:
	for (k=0; k<(dim+1); k++)	{
		free(A[k]);
	}
	free(A);
	free(b);
	free(s);
	*/

	free(ErroVec);
	free(OmegaVec);

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
		free(txxOld[i]);
		free(txyOld[i]);
		free(tyyOld[i]);

		free(dudx[i]);
		free(dudy[i]);
		free(dvdx[i]);
		free(dvdy[i]);

		free(uNew[i]);
		free(vNew[i]);
		free(psiNew[i]);
		free(omegaNew[i]);
		free(pNew[i]);
		free(txxNew[i]);
		free(txyNew[i]);
		free(tyyNew[i]);

		free(VecInd[i]);
	}
	free(uOld);
	free(vOld);
	free(psiOld);
	free(omegaOld);
	free(pOld);
	free(txxOld);
	free(txyOld);
	free(tyyOld);

	free(dudx);
	free(dudy);
	free(dvdx);
	free(dvdy);

	free(uNew);
	free(vNew);
	free(psiNew);
	free(omegaNew);		
	free(pNew);
	free(txxNew);
	free(txyNew);
	free(tyyNew);

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
	txxOld=NULL;
	txyOld=NULL;
	tyyOld=NULL;

	uNew=NULL;
	vNew=NULL;
	psiNew=NULL;
	omegaNew=NULL;
	pNew=NULL;
	txxNew=NULL;
	txyNew=NULL;
	tyyNew=NULL;

	VecInd=NULL;

	return (0);
}
