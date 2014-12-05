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

void QuantumForce::quantumforceOpt (const mat &r , mat &qforce) {
  
  mat jgradient  = zeros(nParticles,nDimensions);

  for(int i = 0; i < nParticles; i++) {
    for(int j = 0; j < nParticles; j++) {
      if(i != j) {
	double rij = norm(r.row(i) - r.row(j));
	jgradient.row(i) += wf->amat(i,j)*(r.row(i) - r.row(j))/rij/((1+beta*rij)*(1+beta*rij));
      }
    }
  }
  
  /*
  mat jgradient = zeros(nParticles, nDimensions);
  for (int k = 0; k < nParticles; k++){
    for (int i = 0; i < k; i++){
      rowvec3 r12Vec = r.row(k) - r.row(i);
      double r12 = sqrt(r12Vec(0)*r12Vec(0) + r12Vec(1)*r12Vec(1) + r12Vec(2)*r12Vec(2));
      jgradient.row(k) += r12Vec*dfdr(r12, k, i)/r12;
    }
    for (int i = k + 1; i < nParticles; i++){
  rowvec3 r12Vec = r.row(i) - r.row(k);
  double r12 = sqrt(r12Vec(0)*r12Vec(0) + r12Vec(1)*r12Vec(1) + r12Vec(2)*r12Vec(2));
  jgradient.row(k) -= r12Vec*dfdr(r12, k, i)/r12;
    }
  }
*/  
  mat sgradient = zeros(nParticles, nDimensions);

  for (int i = 0; i < nParticles/2; i++){
    for (int j = 0; j < nParticles/2; j++){
      sgradient.row(i) += dr->gradient(r,i, j)*wf->inv_up(j, i);
      sgradient.row(i + nParticles/2) += dr->gradient(r,i + nParticles/2, j)*wf->inv_down(j, i);
    }
  }

  qforce = 2*(sgradient + jgradient);
  
  return;

}

double QuantumForce::dfdr(const double &r12, const int &particleNum1, const int &particleNum2){
  return wf->amat(particleNum1, particleNum2)/((1 + beta*r12)*(1 + beta*r12));
}

void QuantumForce::quantumforce (const mat &r , mat &qforce) {

  mat jgradient  = zeros(nParticles,nDimensions);
  rowvec3 g;
  for(int i = 0; i < nParticles; i++) {
    
    for (int n = 0 ; n < nDimensions; n++) {
      g(n) = 0;
    }

    for(int j = 0; j < nParticles; j++) {
      if(i != j) {
	double rij = norm(r.row(i) - r.row(j),2);
	g += wf->amat(i,j)*(r.row(i) - r.row(j))/rij/((1+beta*rij)*(1+beta*rij));
      }
    }
    
    jgradient.row(i) = g;
  }



	double a = alpha;


	double r0 = sqrt(r(0,0)*r(0,0) + r(0,1)*r(0,1) + r(0,2)*r(0,2));
	double r1 = sqrt(r(1,0)*r(1,0) + r(1,1)*r(1,1) + r(1,2)*r(1,2));
	double r2 = sqrt(r(2,0)*r(2,0) + r(2,1)*r(2,1) + r(2,2)*r(2,2));
	double r3 = sqrt(r(3,0)*r(3,0) + r(3,1)*r(3,1) + r(3,2)*r(3,2));



	mat sumat  = zeros( 4,3);

	mat sdmat  = zeros( 4,3);
	sumat[0,0] = a*r(0,0)*((-0.25*a*r0 + 1.0)*exp(a*(r0 + 0.5*r1)) + (0.5*a*r1 - 1.0)*exp(a*(0.5*r0 + r1)))*exp(-1.5*a*(r0 + r1))/r0;
	sdmat[0,0] = 0;
	sumat[0,1] = a*r(0,1)*((-0.25*a*r0 + 1.0)*exp(a*(r0 + 0.5*r1)) + (0.5*a*r1 - 1.0)*exp(a*(0.5*r0 + r1)))*exp(-1.5*a*(r0 + r1))/r0;
	sdmat[0,1] = 0;
	sumat[0,2] = a*r(0,2)*((-0.25*a*r0 + 1.0)*exp(a*(r0 + 0.5*r1)) + (0.5*a*r1 - 1.0)*exp(a*(0.5*r0 + r1)))*exp(-1.5*a*(r0 + r1))/r0;
	sdmat[0,2] = 0;
	sumat[1,0] = a*r(1,0)*(-(0.5*a*r0 - 1.0)*exp(a*(r0 + 0.5*r1)) + (0.25*a*r1 - 1.0)*exp(a*(0.5*r0 + r1)))*exp(-1.5*a*(r0 + r1))/r1;
	sdmat[1,0] = 0;
	sumat[1,1] = a*r(1,1)*(-(0.5*a*r0 - 1.0)*exp(a*(r0 + 0.5*r1)) + (0.25*a*r1 - 1.0)*exp(a*(0.5*r0 + r1)))*exp(-1.5*a*(r0 + r1))/r1;
	sdmat[1,1] = 0;
	sumat[1,2] = a*r(1,2)*(-(0.5*a*r0 - 1.0)*exp(a*(r0 + 0.5*r1)) + (0.25*a*r1 - 1.0)*exp(a*(0.5*r0 + r1)))*exp(-1.5*a*(r0 + r1))/r1;
	sdmat[1,2] = 0;
	sumat[2,0] = 0;
	sdmat[2,0] = a*r(2,0)*((-0.25*a*r2 + 1.0)*exp(a*(r2 + 0.5*r3)) + (0.5*a*r3 - 1.0)*exp(a*(0.5*r2 + r3)))*exp(-1.5*a*(r2 + r3))/r2;
	sumat[2,1] = 0;
	sdmat[2,1] = a*r(2,1)*((-0.25*a*r2 + 1.0)*exp(a*(r2 + 0.5*r3)) + (0.5*a*r3 - 1.0)*exp(a*(0.5*r2 + r3)))*exp(-1.5*a*(r2 + r3))/r2;
	sumat[2,2] = 0;
	sdmat[2,2] = a*r(2,2)*((-0.25*a*r2 + 1.0)*exp(a*(r2 + 0.5*r3)) + (0.5*a*r3 - 1.0)*exp(a*(0.5*r2 + r3)))*exp(-1.5*a*(r2 + r3))/r2;
	sumat[3,0] = 0;
	sdmat[3,0] = a*r(3,0)*(-(0.5*a*r2 - 1.0)*exp(a*(r2 + 0.5*r3)) + (0.25*a*r3 - 1.0)*exp(a*(0.5*r2 + r3)))*exp(-1.5*a*(r2 + r3))/r3;
	sumat[3,1] = 0;
	sdmat[3,1] = a*r(3,1)*(-(0.5*a*r2 - 1.0)*exp(a*(r2 + 0.5*r3)) + (0.25*a*r3 - 1.0)*exp(a*(0.5*r2 + r3)))*exp(-1.5*a*(r2 + r3))/r3;
	sumat[3,2] = 0;
	sdmat[3,2] = a*r(3,2)*(-(0.5*a*r2 - 1.0)*exp(a*(r2 + 0.5*r3)) + (0.25*a*r3 - 1.0)*exp(a*(0.5*r2 + r3)))*exp(-1.5*a*(r2 + r3))/r3;
	double updet = -(-0.5*a*r0 + 1.0)*exp(-0.5*a*r0)*exp(-a*r1) + (-0.5*a*r1 + 1.0)*exp(-a*r0)*exp(-0.5*a*r1);

	double downdet = -(-0.5*a*r2 + 1.0)*exp(-0.5*a*r2)*exp(-a*r3) + (-0.5*a*r3 + 1.0)*exp(-a*r2)*exp(-0.5*a*r3);

	//	for(int j = 0; j < nParticles; j++) {
	//	for(int i = 0; i < nDimensions; i++) {
	// qforce(j,i) = 2*sumat(j,i)/updet + 2*sdmat(j,i)/downdet + 2*gradient(j,i);
	//	}
	//}

	qforce= 2*sumat/updet + 2*sdmat/downdet + 2*jgradient;

	return;
}

