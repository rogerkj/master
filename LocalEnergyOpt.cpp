#include "setup.h"

#include <armadillo>
#include <math.h>
#include "LocalEnergyOpt.h"
#include "Orbitals.h"
using namespace arma;
using namespace std;


LocalEnergyOpt::LocalEnergyOpt(int nDim,int nPart,double _R,int ch,double a,double b,Orbitals* _dr,WaveFunction* _wf) {
	dimention = nDim;
	nParticles = nPart;
	
	R = _R;

	charge = ch;
	alpha = a;
	beta = b;

	dr = _dr;
	wf = _wf;
}

double LocalEnergyOpt::get_local_energy_diatomic(const mat &r) {



  double jlaplacian = 0.0;


  rowvec3 rNuclei; 
  rNuclei << R/2 << 0 << 0;

#ifdef JASTROW

  rowvec3 rijVec;

  mat jgradient  = zeros(nParticles,dimention);
 
  for(int i = 0; i < nParticles; i++) {
    for(int j = 0; j < nParticles; j++) {
      if(i != j) {
	rijVec = r.row(i) - r.row(j);
	double rij = norm(rijVec,2);
	jgradient.row(i) += rijVec*dfdr(rij,i,j)/rij;
	jlaplacian += 2*dfdr(rij,i,j)/rij + d2fdr2(rij,i,j);
      }
    }
  }

  double gradient2 = 0.0;
  for (int i = 0; i < nParticles; i++) {
    gradient2 += dot(jgradient.row(i),jgradient.row(i));
  }
  
  jlaplacian += gradient2;
  
#endif
  
  mat udgradient  = zeros(nParticles,dimention);

  for (int i = 0; i < nParticles/2; i++){
    for (int j = 0; j < nParticles/2; j++){
      udgradient.row(i) += dr->gradient(r, i, j)*wf->inv_up(j, i);
      udgradient.row(i + nParticles/2) += dr->gradient(r, i + nParticles/2, j)*wf->inv_down(j, i);
    }
  }  
  

  double udlaplacian = 0;
  
  for (int i = 0; i < nParticles/2; i++){
    for (int j = 0; j < nParticles/2; j++){
      udlaplacian += dr->laplacian(r,i,j)*wf->inv_up(j,i) + dr->laplacian(r, i + nParticles/2 , j)* wf->inv_down(j,i); 

    }
  }

  double coreEl = 0.0;
  for (int i = 0; i < nParticles; i++) {
    double rn1 = norm(r.row(i) - rNuclei);
    double rn2 = norm(r.row(i) + rNuclei);
    coreEl -= charge*(1/rn1 + 1/rn2);
  }
  
  double electronEl = 0.0;
  for(int i = 0; i < nParticles; i++) {
    for(int j = i+1; j < nParticles; j++) {
      double r12 = norm(r.row(i) - r.row(j));
      electronEl += 1/r12;		}	}
  
  double protonEl = charge*charge/R;


  double gradient = 0.0;
  
#ifdef JASTROW

  for(int i = 0 ; i < nParticles; i++) {
    gradient += dot(udgradient.row(i),jgradient.row(i));
  }
  
#endif

  return -0.5*udlaplacian - 0.5*jlaplacian - gradient + coreEl + electronEl + protonEl;
  
}



double LocalEnergyOpt::get_local_energy(const mat &r) {
   



  double jlaplacian = 0.0;

#ifdef JASTROW

  rowvec3 rijVec;

  mat jgradient  = zeros(nParticles,dimention);
 
  for(int i = 0; i < nParticles; i++) {
    for(int j = 0; j < nParticles; j++) {
      if(i != j) {
	rijVec = r.row(i) - r.row(j);
	double rij = norm(rijVec,2);
	jgradient.row(i) += rijVec*dfdr(rij,i,j)/rij;
	jlaplacian += 2*dfdr(rij,i,j)/rij + d2fdr2(rij,i,j);
      }
    }
  }

  double gradient2 = 0.0;
  for (int i = 0; i < nParticles; i++) {
    gradient2 += dot(jgradient.row(i),jgradient.row(i));
  }
  
  jlaplacian += gradient2;
  
#endif
  
  mat udgradient  = zeros(nParticles,dimention);

  for (int i = 0; i < nParticles/2; i++){
    for (int j = 0; j < nParticles/2; j++){
      udgradient.row(i) += dr->gradient(r, i, j)*wf->inv_up(j, i);
      udgradient.row(i + nParticles/2) += dr->gradient(r, i + nParticles/2, j)*wf->inv_down(j, i);
    }
  }  
  

  double udlaplacian = 0;
  
  for (int i = 0; i < nParticles/2; i++){
    for (int j = 0; j < nParticles/2; j++){
      udlaplacian += dr->laplacian(r,i,j)*wf->inv_up(j,i) + dr->laplacian(r, i + nParticles/2 , j)* wf->inv_down(j,i); 

    }
  }

  double coreEl = 0.0;
  for (int i = 0; i < nParticles; i++) {
    double rn = norm(r.row(i));
    coreEl -= charge/rn;
  }
  
  double electronEl = 0.0;
  for(int i = 0; i < nParticles; i++) {
    for(int j = i+1; j < nParticles; j++) {
      double r12 = norm(r.row(i) - r.row(j));
      electronEl += 1/r12;		}	}
  
  double gradient = 0.0;
  
#ifdef JASTROW

  for(int i = 0 ; i < nParticles; i++) {
    gradient += dot(udgradient.row(i),jgradient.row(i));
  }
  
#endif




  return -0.5*udlaplacian - 0.5*jlaplacian - gradient + coreEl + electronEl;
  
}



double LocalEnergyOpt::dfdr(const double &r12, const int &particleNum1, const int &particleNum2){
  return wf->amat(particleNum1, particleNum2)/((1 + beta*r12)*(1 + beta*r12));
}
double LocalEnergyOpt::d2fdr2(const double &r12, const int &particleNum1, const int &particleNum2){
  return -2*wf->amat(particleNum1, particleNum2)*beta/((1 + beta*r12)*(1 + beta*r12)*(1 + beta*r12));
}
