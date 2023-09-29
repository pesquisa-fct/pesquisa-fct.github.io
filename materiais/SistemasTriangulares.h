#ifndef SISTEMAS_TRIANGULARES_H
#define SISTEMAS_TRIANGULARES_H

void SistemaTriangularInferior(double **A, double *b, double *x, int n, int lu);
void SistemaTriangularSuperior(double **A, double *b, double *x, int n);

#endif
