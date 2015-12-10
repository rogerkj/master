#include "setup.h"

#include <armadillo>
#include <iostream>
#include "QuantumForce.h"
using namespace arma;
using namespace std;
QuantumForce::QuantumForce(int nDim,int nPart,int ch,double a,double b,Orbitals* _dr, WaveFunction* _wf) :
	h(0.001),
	h2(1000000)
{
	nDimensions = nDim;
	nParticles = nPart;
	charge = ch;
	alpha = a;
	beta = b;
	wf = _wf;
	dr = _dr;
}

mat QuantumForce::quantumforceOpt (const mat &r) {
  /*  
  mat jgradient  = zeros(nParticles,nDimensions);

  for(int i = 0; i < nParticles; i++) {
    for(int j = 0; j < nParticles; j++) {
      if(i != j) {
	double rij = norm(r.row(i) - r.row(j));
	jgradient.row(i) += wf->amat(i,j)*(r.row(i) - r.row(j))/rij/((1+beta*rij)*(1+beta*rij));
      }
    }
  }
  */

  //jastrow part
  mat jgradient = zeros(nParticles, nDimensions);
 
#ifdef JASTROW

  for (int k = 0; k < nParticles; k++){
    for (int i = 0; i < k; i++){
      rowvec3 r12Vec = r.row(k) - r.row(i);
      double r12 = sqrt(r12Vec(0)*r12Vec(0) + r12Vec(1)*r12Vec(1) + r12Vec(2)*r12Vec(2));
      
      if (r12 == 0)
	cout << "r12 is zero!" << endl;

      jgradient.row(k) += r12Vec*dfdr(r12, k, i)/r12;
    }
    for (int i = k + 1; i < nParticles; i++){
      rowvec3 r12Vec = r.row(i) - r.row(k);
      double r12 = sqrt(r12Vec(0)*r12Vec(0) + r12Vec(1)*r12Vec(1) + r12Vec(2)*r12Vec(2));

    if (r12 == 0)
	cout << "r12 is zero!" << endl;

      jgradient.row(k) -= r12Vec*dfdr(r12, k, i)/r12;
    }
  }
  

#endif

  //Slater part
  mat sgradient = zeros(nParticles, nDimensions);
  

  for (int i = 0; i < nParticles/2; i++){
    for (int j = 0; j < nParticles/2; j++){
      sgradient.row(i) += dr->gradient(r,i, j)*wf->inv_up(j, i);
      sgradient.row(i + nParticles/2) += dr->gradient(r,i + nParticles/2, j)*wf->inv_down(j, i);
    }
  }

  return 2*(sgradient + jgradient);

}

double QuantumForce::dfdr(double r12, int particleNum1, int particleNum2){
  return wf->amat(particleNum1, particleNum2)/((1 + beta*r12)*(1 + beta*r12));
}

