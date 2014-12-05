#ifndef LOCALENERGY_H
#define LOCALENERGY_H

#include <armadillo>

#include "WaveFunction.h"

using namespace arma;
using namespace std;

class LocalEnergy
{

 public:
  LocalEnergy(int nDim,int nPart,int ch,double a,double b,WaveFunction* _wf );

  double get_local_energy(const mat &r);

 private:

  int dimention;
  int nParticles;

  int charge;

  double alpha;
  double beta;

  WaveFunction* wf;


};

#endif
