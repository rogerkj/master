#ifndef DIATOMIC_H
#define DIATOMIC_H

#include <armadillo>

#include "Orbitals.h"

using namespace arma;
using namespace std;

class Diatomic : public Orbitals {

 public:
  Diatomic(int nDim,int nPart,int ch,double a,double b);

  double waveFunction(const mat &r,int nParticle,int orbital);

  rowvec gradient(const mat &r,int nParticle,int orbital); 

  double laplacian(const mat &r,int nParticle,int orbital);
  
  void setR(double R);

 private:

  mat rNuclei;

  Orbitals* dr;

};

#endif
