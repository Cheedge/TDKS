#include<stdio.h>
#include<math.h>

/* Complex datatype */
struct _dcomplex { double re, im; };
typedef struct _dcomplex dcomplex;

/* ZGEEV prototype */
extern void zgeev( char* jobvl, char* jobvr, int* n, dcomplex* a,
                int* lda, dcomplex* w, dcomplex* vl, int* ldvl, dcomplex* vr, int* ldvr,
                dcomplex* work, int* lwork, double* rwork, int* info );

/*
 * define Crank–Nicholson algorithm
 * e^(-iĤΔτ)≈(1-iĤ∆τ/2)/(1+iĤΔτ/2)
 * (1+i*Ĥ(τ_{j+½})∆τ/2)ψ(τ_{j+1})=(1-i*Ĥ(τ_{j+½})∆τ/2)ψ{τ_j}
 */
int main(int argc, char *argv[])
{
	int i, t;
	double 

/*
 * input1: phi_t=ψ(τ_{j+½}), to construct Ĥ(τ_{j+½})=K(τ_{j+½})+V0(x)+V(τ_{j+½})
 * input2: phi_t0=ψ{τ_j}
 * input3: dt=∆τ
 */
	
	// construct H0 and Ht
	for (i=0; i<M-1; i++)
	{
		K[i][i]=2.0
		K[i][i+1]=-1.0;
		K[i+1][i]=-1.0;
		V0[i][i]=0.5*x[i]*x[i];
//    ft = np.exp(-(t**2)/2)
 //   Vex = constant.A * ft * x * np.sin(constant.omega * t)
		Vt[i][i]=exp(-t*t/2.0) * x[i] * sin(t);
		H[t]=0.5*K/(dx*dx)+V0+Vt;
	}

	// solve eigen problem
	
	//C-N

	return 0;
}
