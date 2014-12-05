#ifndef QUANTUMFORCENUMERIC_H
#define QUANTUMFORCENUMERIC_H

using namespace arma;
using namespace std;

class QuantumForceNumeric
{

 public:

  QuantumForceNumeric(int nDim, int nPart,int ch,double a,double b,WaveFunction* _wf);

  void quantumforce(const mat &r , mat &qforce ,double wfunc );

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





