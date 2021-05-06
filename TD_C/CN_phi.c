/*** THis is 1d CN for Schro Eq ***/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main()
{
 int i,j,k,N,t,NT;
 double *phi,*phit;
 double L,hL,dx,T,dt,x,xp,sum;
 double al,tmp,invdx2; 
 double *D,*E,*EL;
 long int c1,c2,ok; 
  
/*** Define variable ***/

 N=201;

 L= 20.0;
 hL = L/2.0;
 dx = L/(N-1);

 invdx2 = 1.0/(dx*dx);

 T=2;
 NT=2001;
 dt= T/(NT-1);

 al = 1.0;
 
/*** Allocate the memory ***/

 phi = (double *)malloc(sizeof(double)*2*N);
 phit = (double *)malloc(sizeof(double)*2*N);
 

 D = (double *)malloc(sizeof(double)*2*N);
 E = (double *)malloc(sizeof(double)*2*(N-1));
 EL = (double *)malloc(sizeof(double)*2*(N-1));
 
/*** Initialize density matrix  ***/  	 	

 sum = 0.0;

 tmp = pow(al/M_PI,0.25);

 for(i=0;i<N;i++)
  {
   x = i*dx-hL;
   	    
   phi[2*i] = tmp*exp(-0.5*al*x*x);
   phi[2*i+1] = 0.0; 
    
   sum += phi[2*i]*phi[2*i]+phi[2*i+1]*phi[2*i+1]); 
      
  }
	
 sum *= dx;
 
 printf("# Inorm = %3.16f\n",sum);	

/*** Start the Propogation **/

for(t=0;t<NT;t++)
 {
 
 
 for(i=0;i<N;i++)
  {
   phit[2*i] = phi[2*i] -0.25*invdx2*dt*(phi[2*(i+1)+1] - 2.0*phi[2*i+1] + phi[2*(i-1)+1]); 
   phit[2*i+1] = phi[2*i+1] +0.25*invdx2*dt*(phi[2*(i+1)] - 2.0*phi[2*i] + phi[2*(i-1)]);  	
  	
  	  	
   D[2*i] = 1.0;
   D[2*i+1]= 0.5*invdx2*dt;   

  
   if( i < (N-1) )
    {
     E[2*i] = 0.0;
     E[2*i+1] = -0.25*invdx2*dt;
     EL[2*i] = E[2*i];
     EL[2*i+1] = E[2*i+1];
    }
  }
 
 
 c1 = N;
 c2 = 1;

 zgtsv_(&c1,&c2,EL,D,E,phit,&c1,&ok);

 sum = 0.0;

 for(i=0;i<N;i++)
  {
   x = i*dx-hL;
   
   phi[2*i] = phit[2*i];
   phi[2*i+1] = phit[2*i+1];
   
   sum += phi[2*i]*phi[2*i]+phi[2*i+1]*phi[2*i+1]; 
   
/* printf("%lf %3.16f\n",x,phi[2*i]*phi[2*i]+phi[2*i+1]*phi[2*i+1]);*/
   
  }
	
 sum *= dx;
 
 printf("# %lf norm = %3.16f\n",(t+1)*dt,sum);

 }
 

 sum = 0.0;

 for(i=0;i<N;i++)
  {
   x = i*dx-hL;
   
   sum += phi[2*i]*phi[2*i]+phi[2*i+1]*phi[2*i+1]; 
   
   printf("%lf %3.16f\n",x,phi[2*i]*phi[2*i]+phi[2*i+1]*phi[2*i+1]);
   
  }
	
 sum *= dx;
 
 printf("# Fnorm = %3.16f\n",sum); 
 
 
/*** Free the memory ***/
	
 free(phi);	
 free(phit);	
 
 free(D);
 free(E);
 free(EL);

/*** Exit the program ***/

 exit(0);
} 