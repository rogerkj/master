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

  void quantumforce(const mat &r , mat &qforce);

  void quantumforceOpt (const mat &r , mat &qforce);

 private:

  double h,h2;

  int nParticles;
  int nDimensions;

  int charge;

  double alpha;
  double beta;

  WaveFunction* wf;
  Orbitals* dr;

  double dfdr(const double &r12, const int &particleNum1, const int &particleNum2);

};

#endif





