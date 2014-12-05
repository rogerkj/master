#include <armadillo>
#include <math.h>
#include "LocalEnergyOpt.h"
#include "Derivates.h"
using namespace arma;
using namespace std;


LocalEnergyOpt::LocalEnergyOpt(int nDim,int nPart,int ch,double a,double b,Derivates* _dr,WaveFunction* _wf) {
	dimention = nDim;
	nParticles = nPart;
	charge = ch;
	alpha = a;
	beta = b;

	dr = _dr;
	wf = _wf;
}

double LocalEnergyOpt::get_local_energy(const mat &r) {
  
  mat jgradient  = zeros(nParticles,dimention);
    double jlaplacian = 0.0;
  for(int i = 0; i < nParticles; i++) {
    for(int j = 0; j < nParticles; j++) {
      if(i != j) {
	double rij = norm(r.row(i) - r.row(j),2);
	jgradient.row(i) += wf->amat(i,j)*(r.row(i) - r.row(j))/rij/((1+beta*rij)*(1+beta*rij));
	jlaplacian += wf->amat(i,j)*2/(rij*(1+beta*rij)*(1+beta*rij)*(1+beta*rij));
      }
    }
  }
  double gradient2 = 0.0;
  for (int i = 0; i < nParticles; i++) {
    gradient2 += dot(jgradient.row(i),jgradient.row(i));
  }
  
  jlaplacian += gradient2;
  
  mat ugradient  = zeros(nParticles,dimention);
  mat dgradient  = zeros(nParticles,dimention);
  double ulaplacian = 0;
  double dlaplacian = 0;
  
  rowvec3 g,l;
  
  for (int i = 0; i < nParticles ; i++) {
    
    for (int n = 0 ; n < dimention; n++) {
      g(n) = 0;
      l(n) = 0;
    }
  
    for (int j = 0 ; j < nParticles/2; j++) {
      if (i < nParticles/2) {
	g += dr->gradient(r,i,j)*wf->inv_up(j,i);
	l += dr->laplacian(r,i,j)*wf->inv_up(j,i);
      }else {
	g += dr->gradient(r,i,j)*wf->inv_down(j,i-nParticles/2);
	l += dr->laplacian(r,i,j)*wf->inv_down(j,i-nParticles/2);
      }
    }
          
    for (int n = 0 ; n < dimention; n++) {
      if (i < nParticles) {
	ugradient(i,n) = g(n);
	ulaplacian += l(n);
      }else {
	dgradient(i,n) = g(n);
	dlaplacian += l(n);
      }
    }
  }
  
  double coreEl = 0.0;
  for (int i = 0; i < nParticles; i++) {
    double rn = norm(r.row(i));
    coreEl -= charge/rn;
  }
  
  double electronEl = 0.0;
  for(int i = 0; i < nParticles; i++) {
    for(int j = i+1; j < nParticles; j++) {
      double r12 = norm(r.row(i) - r.row(j));
      electronEl += 1/r12;		}	}
  
  double udgradient = 0.0;
  double ujgradient = 0.0;
  double djgradient = 0.0;
  
  for(int i = 0 ; i < nParticles; i++) {
    udgradient += dot(ugradient.row(i),dgradient.row(i));
    ujgradient += dot(jgradient.row(i),ugradient.row(i));
    djgradient += dot(jgradient.row(i),dgradient.row(i));
  }
  
  return -0.5*ulaplacian - 0.5*dlaplacian - 0.5*jlaplacian - udgradient - ujgradient - djgradient + coreEl + electronEl;
  
}
