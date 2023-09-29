#ifndef METODO_POTENCIA_H
#define METODO_POTENCIA_H

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "FuncoesAuxiliares.h"
#include "LUdecomposition.h"
#include "SistemasTriangulares.h"

double MetodoPotencia(double **A, double *u, int n, double tol, int inv);

#endif
