#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<stddef.h>
#include<complex.h>

int main(int argc, char *argv[])
{
	int i, j, N;
	int M;
	double *x, dx;
	_Complex double *psi, *K, *V;
// get M, dx 
	if (argc != 3)
	{
		printf("input as: **.x M, dx\n");
	}
	if (sscanf(argv[1], "%d", &M) != 1)
	{
		printf("give M: number of points\n");
	}
	if (sscanf(argv[2], "%lf", &dx) != 1)
	{
		printf("give dx the delta_x\n");
	}
// allocate variables, arrays.

	psi = (_Complex double*)malloc(sizeof(*psi)*M);/* give M psi*/
	K = (_Complex double*)malloc(sizeof(*K)*M*M);
	V = (_Complex double*)malloc(sizeof(*V)*M*M);

// calculate knetic energy matrix K(x)
//K[i][j] = -(psi[i+1] - 2*psi[i] + psi[i-1])/(2*pow(dx, 2));

	int delta = pow(dx, 2);

	for (i=0; i<=M-1; i++)
	{
		j = i+1
		*(K + i * M + i) = 1.0;
		*(K + i * M + j) = -0.5
		*(K + j + i * M) = -0.5;
	}
	K = K/delta

// calculate potential V(x)
	for (i=0; i<=M-1; i++)
	{
		V[i][i] = pow(x[i], 2)/2.0;
	}

// solve harmitonion


// use Crack-Nilcson to propogate

// solve TDKS

// free arrays.


	return 0;
}
