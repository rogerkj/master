#ifndef ORBITALS_H
#define ORBITALS_H

#include <armadillo>

using namespace arma;
using namespace std;

class Orbitals {

 public:
  Orbitals();

  virtual void setR(double R)=0 ;

  virtual double waveFunction(const mat &r,int nParticle,int orbital)=0;

  virtual rowvec gradient(const mat &r,int nParticle,int orbital)=0; 

  virtual double laplacian(const mat &r,int nParticle,int orbital)=0;
  
  protected:


  int nDimensions;
  int nParticles;

  int charge;

  double alpha;
  double beta;

};

#endif
