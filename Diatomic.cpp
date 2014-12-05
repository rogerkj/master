#include <armadillo>
#include <math.h>

#include "setup.h"
#include "Diatomic.h"
#include "Gaussians.h"
#include "Hydrogenic.h"

using namespace arma;
using namespace std;


Diatomic::Diatomic(int nDim,int nPart,int ch,double a,double b) {
	nDimensions = nDim;
	nParticles = nPart;
	charge = ch;
	alpha = a;
	beta = b;


#ifndef GAUSSIAN

	dr = new Hydrogenic(nDim,nPart,ch,a,b);

#else

	dr = new Gaussians(nDim,nPart,ch,a,b);

#endif

	rNuclei = zeros(nParticles, nDimensions); 

}

void Diatomic::setR(double R){
  for (int i = 0; i < nParticles; i++){
    rNuclei(i,0) = R/2;
  }
}

double Diatomic::waveFunction(const mat &r,int nParticle,int orbital) {

  double value;
 
  if (orbital % 2 == 0){
    value = dr->waveFunction( r + rNuclei, nParticle, orbital/2)
      + dr->waveFunction( r - rNuclei, nParticle, orbital/2);
  } else {
    value = dr->waveFunction( r + rNuclei, nParticle, orbital/2)
      - dr->waveFunction( r - rNuclei, nParticle, orbital/2);
  }
  return value;

}

rowvec Diatomic::gradient(const mat &r,int nParticle,int orbital) {

  rowvec gradient(nDimensions);

  if (orbital % 2 == 0){
    gradient = dr->gradient( r + rNuclei, nParticle, orbital/2)
      + dr->gradient( r - rNuclei, nParticle, orbital/2);
  } else {
    gradient = dr->gradient( r + rNuclei, nParticle, orbital/2)
      - dr->gradient( r - rNuclei, nParticle, orbital/2);
  }
  return gradient;
}
double Diatomic::laplacian(const mat &r,int nParticle,int orbital) {

  double laplacian;

  if (orbital % 2 == 0){
    laplacian = dr->laplacian( r + rNuclei, nParticle, orbital/2)
      + dr->laplacian( r - rNuclei,nParticle, orbital/2);
  } else {
    laplacian = dr->laplacian( r + rNuclei, nParticle, orbital/2)
      - dr->laplacian( r - rNuclei, nParticle, orbital/2);
  }
  
  return laplacian;
  
}
