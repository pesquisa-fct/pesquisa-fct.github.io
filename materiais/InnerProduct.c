#include "InnerProduct.h"

double InnerProduct(double *x, double *y, int dim){

	int i;
	double sum=0.0;
	for (i=0; i<dim; i++){
		sum+= x[i]*y[i];
	}

	return (sum);
}
