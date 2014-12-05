#ifndef LOCALENERGYNUMERIC_H
#define LOCALENERGYNUMERIC_H

#include <armadillo>

#include "WaveFunction.h"

using namespace arma;
using namespace std;

class LocalEnergyNumeric
{

 public:

  LocalEnergyNumeric(int nDim, int nPart,int ch,double a,double b,WaveFunction* _wf);

  double get_local_energy(const mat &r);

 private:

  double h,h2;

  WaveFunction* wf;

  int nParticles;
  int nDimensions;

  int charge;

  double alpha;
  double beta;


};

#endif
