#ifndef LOCALENERGYOPT_H
#define LOCALENERGYOPT_H

#include <armadillo>

#include "LocalEnergyOpt.h"
#include "Orbitals.h"
#include "WaveFunction.h"

using namespace arma;
using namespace std;

class LocalEnergyOpt
{

 public:
  LocalEnergyOpt(int nDim,int nPart,double _R,int ch,double a,double b,Orbitals* _dr,WaveFunction* _wf);

  double get_local_energy(const mat &r);

  double get_local_energy_diatomic(const mat &r);

 private:

  int dimention;
  int nParticles;

  double R;

  int charge;

  double alpha;
  double beta;

  Orbitals* dr;
  WaveFunction* wf;

  double dfdr(const double &r12, const int &particleNum1, const int &particleNum2);
 
  double d2fdr2(const double &r12, const int &particleNum1, const int &particleNum2);

};

#endif
