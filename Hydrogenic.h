#ifndef HYDROGENIC_H
#define HYDROGENIC_H

#include <armadillo>

#include "Orbitals.h"

using namespace arma;
using namespace std;

class Hydrogenic : public Orbitals 
{

 public:
  Hydrogenic(int nDim,int nPart,int ch,double a,double b);

  void setR(double R);

  double waveFunction(const mat &r,int nParticle,int orbital);

  rowvec gradient(const mat &r,int nParticle,int orbital); 

  double laplacian(const mat &r,int nParticle,int orbital);
  

};

#endif
