#ifndef FUNCOES_AUXILIARES_H
#define FUNCOES_AUXILIARES_H

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double Norm(double *v, int n);
double MaxNorm(double *v, int n);
double InnerProduct(double *u, double *v, int n);
double *MatrixProduct(double **A, double *v, int n);
void MatrixProductPenta(double Centro, double DirEsq, double CimaBaixo, double *v, double *prod, int Ny, int n);
void MenorMatriz(double **A, double **M, int i, int j, int n);
double Determinant(double **A, int n);
void MatrizDeslocada(double **A, double **Aq, int n, double q);


#endif
