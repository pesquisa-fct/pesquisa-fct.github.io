#ifndef CONJUGATE_GRADIENT_METHOD_H
#define CONJUGATE_GRADIENT_METHOD_H

#include<stdio.h>
#include<stdlib.h>
#include "FuncoesAuxiliares.h"

int ConjugateGradientMethod(double **A, double *b, double *x, int n, double tol);

#endif
