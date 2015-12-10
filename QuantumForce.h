#ifndef QUANTUMFORCE_H
#define QUANTUMFORCE_H

#include "WaveFunction.h"
#include "Orbitals.h"

using namespace arma;
using namespace std;

class QuantumForce
{

 public:

  QuantumForce(int nDim, int nPart,int ch,double a,double b,Orbitals* _dr,WaveFunction* _wf);

  mat quantumforceOpt (const mat &r);

 private:

  double h,h2;

  int nParticles;
  int nDimensions;

  int charge;

  double alpha;
  double beta;

  WaveFunction* wf;
  Orbitals* dr;

  double dfdr(double r12,int particleNum1, int particleNum2);

};

#endif





