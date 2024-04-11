#ifndef VISUALIZACAO_H
#define VISUALIZACAO_H

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

void ImprimeArquivoVTK(double **u, double **v, double **p, double **psi, double **omega, double x0, double xf, double y0, double yf, double *x, double *y, int Nx, int Ny, int n);

#endif
