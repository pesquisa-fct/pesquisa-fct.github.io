#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void **AlocaMatriz(int Linhas, int Colunas, size_t TamanhoElemento, void *ValorInicial);
void LerArquivoVTK(const char *NomeArquivo, double ***u, double ***v, double ***p, double ***psi, double ***omega, int *Nx, int *Ny, double **x, double **y, double *dx, double *dy);
void asymptotics(double **u, double **v, double **p, double **psi, double **omega, int Nx, int Ny, double *x, double *y, int Geometria, int Malha);

int main()
{
    	int Nx, Ny, i, j;
    	double **u, **v, **p, **psi, **omega, *x, *y;    
	double dx, dy;
    	int Geometria, Malha;
    	double Re;
    	Re=1.0;
	//_______________
	//Geometria:
	//_______________
	//4Lx4L:	0
	//8Lx8L:	1
	//16Lx16L:	2
	Geometria = 0;	
	//_______________________
	//Malha:
	//_______________________
	//M1 (dx=dy=0.1):	0
	//M2 (dx=dy=0.05):	1
	//M3 (dx=dy=0.025):	2
	Malha = 0;

	//printf("Antes de leitura\n");
    	LerArquivoVTK("../Results/4Lx4L_M2_Re1/Dados-N10000.vtk", &u, &v, &p, &psi, &omega, &Nx, &Ny, &x, &y, &dx, &dy);    	    

	x[0] = 0.0;
	y[0] = 0.0;
	for (i=1; i<=Nx; i++){
		x[i] = x[i-1] + dx;
	}
	for (j=1; j<=Ny; j++){
		y[j] = y[j-1] + dy;
	}

	/*
	for (i=0; i<=Nx; i++){
		for (j=0; j<=Ny; j++){
			printf("i=%d, j=%d\n", i, j);
			printf("u=%f, v=%f, p=%f, psi=%f, omega=%f\n", u[i][j], v[i][j], p[i][j], psi[i][j], omega[i][j]);getchar();
		}
	}
	*/

	//printf("Apos leitura\n");
    	asymptotics(u, v, p, psi, omega, Nx, Ny, x, y, Geometria, Malha);   
    	
	return 0;
}

void LerArquivoVTK(const char *NomeArquivo, double ***u, double ***v, double ***p, double ***psi, double ***omega, int *Nx, int *Ny, double **x, double **y, double *dx, double *dy)
{
    FILE *arq;
    int i, j, temp;
    double valorInicial;
    double uBaixo, uCima, vEsq, vDir;

    //Criando os arquivos
    arq = fopen(NomeArquivo, "rt");

    if( arq==NULL ) {
        printf("\nProblema na abertura do arquivo.\n");
        exit(0);
    }

    fscanf(arq, "# vtk DataFile Version 2.0\nCavidade\nASCII\n");
    fscanf(arq, "DATASET STRUCTURED_POINTS\n");
    fscanf(arq, "\n");
    fscanf(arq, "DIMENSIONS %d %d 1\n", Nx, Ny);
    (*Nx)--;
    (*Ny)--;
    fscanf(arq, "ORIGIN 0.0 0.0 0.0\n");
    fscanf(arq, "SPACING %lf %lf 1.0\n", dx, dy);
    //fscanf(arq, "SPACING 0.1 0.1 1.0\n");
    //printf("dx=%f, dy=%f\n", *dx, *dy);
    //printf("after dimensinos: Nx=%d, Ny=%d\n", (*Nx), (*Ny));getchar();
    fscanf(arq, "\n");

    //Alocando memoria
    

    valorInicial = 0.0;
    *u = (double **)AlocaMatriz((*Nx)+1, (*Ny)+1, sizeof(double), &valorInicial);
    *v = (double **)AlocaMatriz((*Nx)+1, (*Ny)+1, sizeof(double), &valorInicial);
    *p = (double **)AlocaMatriz((*Nx)+1, (*Ny)+1, sizeof(double), &valorInicial);
    *psi = (double **)AlocaMatriz((*Nx)+1, (*Ny)+1, sizeof(double), &valorInicial);
    *omega = (double **)AlocaMatriz((*Nx)+1, (*Ny)+1, sizeof(double), &valorInicial);
    *x = (double *)malloc( ((*Nx)+1)*sizeof(double) );
    *y = (double *)malloc( ((*Ny)+1)*sizeof(double) );

     //Velocity vectors
    fscanf(arq, "POINT_DATA %d\n", &temp);
    //printf("temp=%d\n", temp);getchar();

    fscanf(arq, "VECTORS velocidade float\n");
    for( j=0; j<=*Ny; j++ ){
        for( i=0; i<=*Nx; i++ ) {
            fscanf(arq, "%lf %lf 0.0000\n", &(*u)[i][j], &(*v)[i][j]);
        }
    }
    fscanf(arq, "\n");

    fscanf(arq, "SCALARS velocidade-u float 1\n");
    fscanf(arq, "LOOKUP_TABLE default\n");
    for(j=0; j<=*Ny; j++) {
        for (i=0; i<=*Nx; i++) {
                fscanf(arq, "%lf\n", &(*u)[i][j]);
		//printf("(*u)[%d][%d]=%f\n", i, j, (*u)[i][j]);//getchar();

        }
    }
    fscanf(arq, "\n");

    fscanf(arq, "SCALARS velocidade-v float 1\n");
    fscanf(arq, "LOOKUP_TABLE default\n");
    for(j=0; j<=*Ny; j++) {
        for(i=0; i<=*Nx; i++) {
                fscanf(arq, "%lf\n", &(*v)[i][j]);
        }
    }
    fscanf(arq, "\n");

    fscanf(arq, "SCALARS vorticity float 1\n");
    fscanf(arq, "LOOKUP_TABLE default\n");
    for(j=0; j<=*Ny; j++)
        for(i=0; i<=*Nx; i++)
            fscanf(arq, "%lf\n", &(*omega)[i][j]);
    fscanf(arq, "\n");

    fscanf(arq, "SCALARS streamfunction float 1\n");
    fscanf(arq, "LOOKUP_TABLE default\n");
    for(j=0; j<=*Ny; j++)
        for(i=0; i<=*Nx; i++)
            fscanf(arq, "%lf\n", &(*psi)[i][j]);
    fscanf(arq, "\n");

    fscanf(arq, "SCALARS pressao float 1\n");
    fscanf(arq, "LOOKUP_TABLE default\n");
    for(j=0; j<=*Ny; j++)
        for(i=0; i<=*Nx; i++)
            fscanf(arq, "%lf\n", &(*p)[i][j]);
    fscanf(arq, "\n");

    fclose(arq);

    return;
}

void **AlocaMatriz(int Linhas, int Colunas, size_t TamanhoElemento, void *ValorInicial)
{
    void **matriz;
    int i, j, byteAtualValor, byteAtualMatriz;
    unsigned char *charValorInicial, *charMatriz;

    charValorInicial = (unsigned char *)ValorInicial;

    matriz = (void **)malloc( Linhas*sizeof(void *) );
    for( i=0; i<Linhas; i++ )
        matriz[i] = malloc( Colunas*TamanhoElemento );

    //Inicializando todos os elementos com o valor inicial (soh funciona se for um tipo int)
    if( ValorInicial!=NULL ) {
        for(i=0; i<Linhas; i++) {
            charMatriz = (unsigned char *)matriz[i];
            byteAtualMatriz = 0;
            for(j=0; j<Colunas; j++) {
                for( byteAtualValor = 0; byteAtualValor<TamanhoElemento; byteAtualValor++ )
                    charMatriz[byteAtualMatriz++] = charValorInicial[byteAtualValor];
            }
        }
    }


    return matriz;
}

void asymptotics(double **u, double **v, double **p, double **psi, double **omega, int Nx, int Ny, double *x, double *y, int Geometria, int Malha)
{
    	int i, j, im, jmDown;
        FILE *pt;
    	double r, xm, ymDown;

    	switch(Geometria){
		case 0://4Lx4L
 			xm=4.0;
			ymDown=1.5;
	    	break;
	    	case 1://8Lx8L
 			xm=8.0;
			ymDown=3.0;
	    	break;
	    	case 2://16Lx16L
 			xm=16.0;
			ymDown=3.0;
	    	break;
       	}

    	//Taking the value of i - middle of channel in x direction
    	for (i=0; i<=Nx; i++)
    	{
		//printf("x[%d]=%f\n", i, x[i]);
        	if ( (fabs(x[i]) - xm) < 1.0e-06 )
        	{
            		im = i;
        	}
    	}
    
    	//Taking the values of j - middle of channel in y direction
    	for (j=0; j<=Ny; j++)
    	{
        	if ( (fabs(y[j]) - ymDown) < 1.0e-06 )
        	{
            		jmDown = j;
        	}
    	}

	pt = fopen("../Results/4Lx4L_M2_Re1/asymptotics.dat", "w");

    	j=jmDown;
    	i=im;
    	r = fabs(y[j] - ymDown);

	//printf("jmDown=%d, im=%d\n", jmDown, im);
	//printf("xm=%f, ym=%f\n", x[i], y[j]);getchar();
    	while ( r <= 1.0 )
    	{
		//printf("r=%f\n", r);
        	//printf("%lf \t %lf \t %lf \t %lf \t %lf \t %lf \n", r, u[i][j], v[i][j], p[i][j], psi[i][j], omega[i][j]);getchar();

        	fprintf(pt, "%lf \t %lf \t %lf \t %lf \t %lf \t %lf \n", r, u[i][j], v[i][j], p[i][j], psi[i][j], omega[i][j]);
		//printf("after\n");
        	j++;
        	r = fabs(y[j] - ymDown);
    	}

    	fclose(pt);

    	//r:	1
    	//u:	2
    	//v:	3
    	//p:	4
    	//psi:	5
    	//omega:6

}

