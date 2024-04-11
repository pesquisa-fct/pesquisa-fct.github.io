#include "Visualizacao.h"

void ImprimeArquivoVTK(double **u, double **v, double **p, double **psi, double **omega, double x0, double xf, double y0, double yf, double *x, double *y, int Nx, int Ny, int n){

    FILE *arq;
    int i, j;
    double dx, dy;
    char nomeArq[600];

    dx = (xf-x0) / Nx;
    dy = (yf-y0) / Ny;

    //Criando os arquivos
    sprintf(nomeArq, "../Results/Dados-N%d.vtk", n);

    //open output file
    arq = fopen(nomeArq, "wt");
    if( !arq ){
        printf("\n\n Problema no fopen do arquivo. \n\n");
        exit(1);
    }

    //vtk type file
    printf("Writing vtk fil...\n");
    fprintf(arq, "# vtk DataFile Version 2.0\nCavidade\nASCII\n");
    fprintf(arq, "DATASET STRUCTURED_POINTS\n");
    fprintf(arq, "\n");
    fprintf(arq, "DIMENSIONS %d %d 1\n", Nx+1, Ny+1);
    fprintf(arq, "ORIGIN 0.0 0.0 0.0\n");
    fprintf(arq, "SPACING %lf %lf 1.0\n", dx, dy);
    fprintf(arq, "\n");

    //Velocity vectors
    fprintf(arq, "POINT_DATA %d\n", (Nx+1)*(Ny+1));
    fprintf(arq, "VECTORS velocidade float\n");
    for( j=0; j<=Ny; j++ ){
        for( i=0; i<=Nx; i++ ) {
            fprintf(arq, "%e %e 0.0000\n", u[i][j], v[i][j] );
        }
    }
    fprintf(arq, "\n");

    fprintf(arq, "SCALARS velocidade-u float 1\n");
    fprintf(arq, "LOOKUP_TABLE default\n");
    for(j=0; j<=Ny; j++) {
        for (i=0; i<=Nx; i++) {
                fprintf(arq, "%e\n", u[i][j]);
        }
    }
    fprintf(arq, "\n");

     fprintf(arq, "SCALARS velocidade-v float 1\n");
    fprintf(arq, "LOOKUP_TABLE default\n");
    for(j=0; j<=Ny; j++) {
        for(i=0; i<=Nx; i++) {
                fprintf(arq, "%e\n", v[i][j]);
        }
    }
    fprintf(arq, "\n");

    fprintf(arq, "SCALARS vorticity float 1\n");
    fprintf(arq, "LOOKUP_TABLE default\n");
    for(j=0; j<=Ny; j++)
        for(i=0; i<=Nx; i++)
            fprintf(arq, "%e\n", omega[i][j]);
    fprintf(arq, "\n");

    fprintf(arq, "SCALARS streamfunction float 1\n");
    fprintf(arq, "LOOKUP_TABLE default\n");
    for(j=0; j<=Ny; j++)
        for(i=0; i<=Nx; i++)
            fprintf(arq, "%e\n", psi[i][j]);
    fprintf(arq, "\n");

    fprintf(arq, "SCALARS pressao float 1\n");
    fprintf(arq, "LOOKUP_TABLE default\n");
    for(j=0; j<=Ny; j++)
        for(i=0; i<=Nx; i++)
            fprintf(arq, "%e\n", p[i][j]);
    fprintf(arq, "\n");

    fclose(arq);

    return;
}

