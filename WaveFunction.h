#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <armadillo>

#include "Orbitals.h"

using namespace arma;


class WaveFunction
{

 public:
  
  WaveFunction(int nPart,int nDim,double a,double b,Orbitals* _dr);
 
  double amat(int i, int j);

  double slater(const mat &r);

  double jastrowOpt(const mat &rNew, const mat &rOld,int i);
  double slaterOpt(const mat &rNew,const mat &rOld,int i);

  double jastrowbrute(const mat &r_new,const mat &r_old);
  double slaterbrute(const mat &rNew,const mat &rOld);


  mat slater_up(const mat &r);
  mat slater_down(const mat &r);

  double waveFunction(const mat &r,int nParticle,int orbital);

  void update_inverse(double Rs,const mat &r,int i);
  void set_inverse(const mat &r);

  mat inv_up;
  mat inv_down;

 private:

  Orbitals* dr;

  double alpha;
  double beta;

  int nDimensions;
  int nParticles;

  double jastrowFactor(const mat &r,int i, int j);

  bool* spin;

};

#endif
