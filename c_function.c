#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <math.h>

double min(double x,double y);
SEXP compute_integral(SEXP phi_G,SEXP fI1,SEXP fI2,SEXP Delta);
double compute_integral_internal(double* phi_G,double* fI1,double* fI2,double Delta,int N,int Ninner);
double* laguerre_density_pos(double* z,int n,double* theta,int m);
double inc_likelihoodC(double x1,double x2,double w,double rC,double* theta1,double* theta2,int m1,int m2,int N);
SEXP whole_inc_likelihood(SEXP Rtheta1,SEXP Rtheta2,SEXP x1,SEXP x2,SEXP w,SEXP rC,SEXP N);


SEXP compute_integral(SEXP phi_G,SEXP fI1,SEXP fI2,SEXP Delta)
{
  // Variables
  double out;
  int i,j;
  int N,Ninner;
  SEXP integral;
  
  N=LENGTH(fI2);
  Ninner=LENGTH(fI1);
  
  // Compute Integral
  out=0.0;
  for(j=0;j<=Ninner-1;j++) {
    for(i=0;i+j<=N-1;i++) {
      out=out+asReal(Delta)*asReal(Delta)*REAL(phi_G)[i]*REAL(fI1)[j]*REAL(fI2)[i+j];
    }
  }
  
  // Create Output
  integral=PROTECT(allocVector(REALSXP,1));
  REAL(integral)[0]=out;
  UNPROTECT(1);
  
  return(integral);
}

SEXP whole_inc_likelihood(SEXP Rtheta1,SEXP Rtheta2,SEXP x1,SEXP x2,SEXP w,SEXP rC,SEXP RN)
{
  int i,n,m1,m2,N;
  double logL;
  double* theta1;
  double* theta2;
  double step;
  SEXP out;
  
  n=LENGTH(x1);
  m1=LENGTH(Rtheta1)-1;
  m2=LENGTH(Rtheta2)-1;
  N=INTEGER(RN)[0];
  
  // Copy thetas to doubles
  theta1=malloc(sizeof(double)*(m1+1));
  theta2=malloc(sizeof(double)*(m2+1));
  for(i=0;i<=m1;i++)
    theta1[i]=REAL(Rtheta1)[i];
  for(i=0;i<=m2;i++)
    theta2[i]=REAL(Rtheta2)[i];
  
  // Compute Likelihood
  logL=0;
  for(i=0;i<=n-1;i++)
  {
    step=inc_likelihoodC(REAL(x1)[i],REAL(x2)[i],REAL(w)[i],REAL(rC)[i],theta1,theta2,m1,m2,N);
    logL=logL+step;
  }
    

  // Create output
  out=PROTECT(allocVector(REALSXP,1));
  REAL(out)[0]=logL;
  UNPROTECT(1);
  
  return(out);
}

double compute_integral_internal(double* phi_G,double* fI1,double* fI2,double Delta,int N,int Ninner)
{
  // Variables
  double out;
  int i,j;

//  N=LENGTH(fI2);
//  Ninner=LENGTH(fI1);
  
  // Compute Integral
  out=0.0;
  for(j=0;j<=Ninner-1;j++)
  {
    for(i=0;i+j<=N-1;i++)
    {
      out=out+Delta*Delta*phi_G[i]*fI1[j]*fI2[i+j];
    }
  }
  
  return(out);
}

// m1=length(theta1)-1, similarly m2
double inc_likelihoodC(double x1,double x2,double w,double rC,double* theta1,double* theta2,int m1,int m2,int N)
{
  double* y;
  double Delta;
  int lim1;
  int i;
  double* ax1;
  double* ax2;
  double* fIx1;
  double* fIx2;
  double* fG;
  double integral;
  
  // Allocate memory
  Delta=x2/((double)(N-1));
  lim1=(int)floor(min(x1,w)/Delta);
  
  y=malloc(sizeof(double)*N);
  ax1=malloc(sizeof(double)*(lim1+1));
  ax2=malloc(sizeof(double)*N);
  
  
  // Approximating sequence for y
  for(i=0;i<=N-1;i++)
    y[i]=Delta*(double)i;
  
  // Approximating sequence for x1 and x2
  for(i=0;i<=lim1;i++)
    ax1[i]=x1-(double)i*Delta;
  
  for(i=0;i<=N-1;i++)
    ax2[i]=x2-(double)i*Delta;

  // Compute Laguerre polynomials on grid
  fIx1=laguerre_density_pos(ax1,lim1+1,theta1,m1);
  for(i=0;i<=lim1;i++)  
    fIx1[i]=fIx1[i]*exp(-rC*(w-x1+ax1[i]));
  
  fIx2=laguerre_density_pos(ax2,N,theta1,m1);
  fG=laguerre_density_pos(y,N,theta2,m2);

  
  // Compute integral
  integral=compute_integral_internal(fG,fIx1,fIx2,Delta,N,lim1+1);

  // Free memory
  free(y);
  free(ax1);
  free(ax2);
  free(fIx1);
  free(fIx2);
  free(fG);
        
  return(log(integral));
}

// n is length of z and length of theta minus 1=m
double* laguerre_density_pos(double* z,int n,double* theta,int m)
{
  double* f;
  double* B;
  double* Zp;
  double* Lp;
  int k,i,j;
  
  f=malloc(sizeof(double)*n);
  Lp=malloc(sizeof(double)*n);
  B=malloc(sizeof(double)*(m+1)*(m+1));
  Zp=malloc(sizeof(double)*n*(m+1));
  
  // Compute the matrix of binomial coefficients
  for(k=1;k<=m+1;k++)
  {
    for(i=0;i<=m;i++)
    {
      if(i==0)
        B[(k-1)+i*(m+1)]=1.0;
      else if(i<=k-1)
        B[(k-1)+i*(m+1)]=(double)(k-1-i+1)/((double)i)*B[(k-1)+(i-1)*(m+1)];
      else
        B[(k-1)+i*(m+1)]=0;
    }
  }
  
  // Compute the powers of z
  for(k=0;k<=n-1;k++)
    Zp[k]=1.0;
  
  if(m+1>=2)
  {
    for(i=1;i<=m;i++)
      for(k=0;k<=n-1;k++)
          Zp[k+i*n]=Zp[k+(i-1)*n]*(-z[k])/((double)i);
  }
      
  // Compute Laguerre Polynomials
  for(k=0;k<=n-1;k++)
  {
    Lp[k]=0;
    for(i=0;i<=m;i++)
      for(j=0;j<=m;j++)
        Lp[k]=Lp[k]+Zp[k+j*n]*B[i+j*(m+1)]*theta[i];
  }
        
  // Compute the density
  for(k=0;k<=n-1;k++)
    f[k]=exp(-z[k])*Lp[k]*Lp[k];
  
  free(Lp);
  free(B);
  free(Zp);
        
  return(f);
}

double min(double x,double y)
{
  if(x<=y)
    return(x);
  else
    return(y);
}