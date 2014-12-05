#include <armadillo>
#include <math.h>
#include "LocalEnergy.h"
using namespace arma;
using namespace std;


LocalEnergy::LocalEnergy(int nDim,int nPart,int ch,double a,double b,WaveFunction* _wf ) {
	dimention = nDim;
	nParticles = nPart;
	charge = ch;
	alpha = a;
	beta = b;

	wf = _wf;

}



double LocalEnergy::get_local_energy(const mat &r) {

	double a = alpha;

	mat gradient  = zeros(nParticles,dimention);
	rowvec g(dimention);

	double laplacian = 0.0;
	for(int i = 0; i < nParticles; i++) {
	  g(0) = 0;
	  g(1) = 0;
	  g(2) = 0;

	  for(int j = 0; j < nParticles; j++) { 
	    if(i != j) {
	      double rij = norm(r.row(i) - r.row(j),2);
	      g += wf->amat(i,j)*(r.row(i) - r.row(j))/rij/((1+beta*rij)*(1+beta*rij));

	      laplacian += wf->amat(i,j)*2/(rij*(1+beta*rij)*(1+beta*rij)*(1+beta*rij));
	     
	    }
	  }

	  gradient(i,0) = g(0);
	  gradient(i,1) = g(1);
	  gradient(i,2) = g(2);
	}
	
	double gradient2 = 0.0;
	
	for (int i = 0; i < nParticles; i++) {
	  gradient2 += dot(gradient.row(i),gradient.row(i));
	}

	laplacian += gradient2;

	double r0 = sqrt(r(0,0)*r(0,0) + r(0,1)*r(0,1) + r(0,2)*r(0,2));
	double r1 = sqrt(r(1,0)*r(1,0) + r(1,1)*r(1,1) + r(1,2)*r(1,2));
	double r2 = sqrt(r(2,0)*r(2,0) + r(2,1)*r(2,1) + r(2,2)*r(2,2));
	double r3 = sqrt(r(3,0)*r(3,0) + r(3,1)*r(3,1) + r(3,2)*r(3,2));
	double d2u = a*(a*r0*r0*r0*r1*(0.5*a*r0 - 1.0)*(r(1,0)*r(1,0) + r(1,1)*r(1,1) + r(1,2)*r(1,2))*exp(a*(8.5*r0 + 8.0*r1)) + a*r0*r0*r0*r1*(r(1,0)*r(1,0)*(-0.125*a*r1 + 0.25) + 0.5*r(1,0)*r(1,0) + r(1,1)*r(1,1)*(-0.125*a*r1 + 0.25) + 0.5*r(1,1)*r(1,1) + r(1,2)*r(1,2)*(-0.125*a*r1 + 0.25) + 0.5*r(1,2)*r(1,2))*exp(a*(8.0*r0 + 8.5*r1)) + a*r0*r1*r1*r1*(-0.5*a*r1 + 1.0)*(r(0,0)*r(0,0) + r(0,1)*r(0,1) + r(0,2)*r(0,2))*exp(a*(8.0*r0 + 8.5*r1)) + a*r0*r1*r1*r1*(r(0,0)*r(0,0)*(0.125*a*r0 - 0.25) - 0.5*r(0,0)*r(0,0) + r(0,1)*r(0,1)*(0.125*a*r0 - 0.25) - 0.5*r(0,1)*r(0,1) + r(0,2)*r(0,2)*(0.125*a*r0 - 0.25) - 0.5*r(0,2)*r(0,2))*exp(a*(8.5*r0 + 8.0*r1)) + r0*r0*r0*r1*r1*(-1.5*a*r0 + 3.0)*exp(a*(8.5*r0 + 8.0*r1)) + r0*r0*r0*r1*r1*(0.75*a*r1 - 3.0)*exp(a*(8.0*r0 + 8.5*r1)) + r0*r0*r0*(0.5*a*r0 - 1.0)*(r(1,0)*r(1,0) + r(1,1)*r(1,1) + r(1,2)*r(1,2))*exp(a*(8.5*r0 + 8.0*r1)) + r0*r0*r0*(r(1,0)*r(1,0)*(-0.25*a*r1 + 0.5) + 0.5*r(1,0)*r(1,0) + r(1,1)*r(1,1)*(-0.25*a*r1 + 0.5) + 0.5*r(1,1)*r(1,1) + r(1,2)*r(1,2)*(-0.25*a*r1 + 0.5) + 0.5*r(1,2)*r(1,2))*exp(a*(8.0*r0 + 8.5*r1)) + r0*r0*r1*r1*r1*(-0.75*a*r0 + 3.0)*exp(a*(8.5*r0 + 8.0*r1)) + r0*r0*r1*r1*r1*(1.5*a*r1 - 3.0)*exp(a*(8.0*r0 + 8.5*r1)) - r1*r1*r1*(0.5*a*r1 - 1.0)*(r(0,0)*r(0,0) + r(0,1)*r(0,1) + r(0,2)*r(0,2))*exp(a*(8.0*r0 + 8.5*r1)) + r1*r1*r1*(r(0,0)*r(0,0)*(0.25*a*r0 - 0.5) - 0.5*r(0,0)*r(0,0) + r(0,1)*r(0,1)*(0.25*a*r0 - 0.5) - 0.5*r(0,1)*r(0,1) + r(0,2)*r(0,2)*(0.25*a*r0 - 0.5) - 0.5*r(0,2)*r(0,2))*exp(a*(8.5*r0 + 8.0*r1)))*exp(-9.0*a*(r0 + r1))/(r0*r0*r0*r1*r1*r1);

	double d2d = a*(a*r2*r2*r2*r3*(0.5*a*r2 - 1.0)*(r(3,0)*r(3,0) + r(3,1)*r(3,1) + r(3,2)*r(3,2))*exp(a*(8.5*r2 + 8.0*r3)) + a*r2*r2*r2*r3*(r(3,0)*r(3,0)*(-0.125*a*r3 + 0.25) + 0.5*r(3,0)*r(3,0) + r(3,1)*r(3,1)*(-0.125*a*r3 + 0.25) + 0.5*r(3,1)*r(3,1) + r(3,2)*r(3,2)*(-0.125*a*r3 + 0.25) + 0.5*r(3,2)*r(3,2))*exp(a*(8.0*r2 + 8.5*r3)) + a*r2*r3*r3*r3*(-0.5*a*r3 + 1.0)*(r(2,0)*r(2,0) + r(2,1)*r(2,1) + r(2,2)*r(2,2))*exp(a*(8.0*r2 + 8.5*r3)) + a*r2*r3*r3*r3*(r(2,0)*r(2,0)*(0.125*a*r2 - 0.25) - 0.5*r(2,0)*r(2,0) + r(2,1)*r(2,1)*(0.125*a*r2 - 0.25) - 0.5*r(2,1)*r(2,1) + r(2,2)*r(2,2)*(0.125*a*r2 - 0.25) - 0.5*r(2,2)*r(2,2))*exp(a*(8.5*r2 + 8.0*r3)) + r2*r2*r2*r3*r3*(-1.5*a*r2 + 3.0)*exp(a*(8.5*r2 + 8.0*r3)) + r2*r2*r2*r3*r3*(0.75*a*r3 - 3.0)*exp(a*(8.0*r2 + 8.5*r3)) + r2*r2*r2*(0.5*a*r2 - 1.0)*(r(3,0)*r(3,0) + r(3,1)*r(3,1) + r(3,2)*r(3,2))*exp(a*(8.5*r2 + 8.0*r3)) + r2*r2*r2*(r(3,0)*r(3,0)*(-0.25*a*r3 + 0.5) + 0.5*r(3,0)*r(3,0) + r(3,1)*r(3,1)*(-0.25*a*r3 + 0.5) + 0.5*r(3,1)*r(3,1) + r(3,2)*r(3,2)*(-0.25*a*r3 + 0.5) + 0.5*r(3,2)*r(3,2))*exp(a*(8.0*r2 + 8.5*r3)) + r2*r2*r3*r3*r3*(-0.75*a*r2 + 3.0)*exp(a*(8.5*r2 + 8.0*r3)) + r2*r2*r3*r3*r3*(1.5*a*r3 - 3.0)*exp(a*(8.0*r2 + 8.5*r3)) - r3*r3*r3*(0.5*a*r3 - 1.0)*(r(2,0)*r(2,0) + r(2,1)*r(2,1) + r(2,2)*r(2,2))*exp(a*(8.0*r2 + 8.5*r3)) + r3*r3*r3*(r(2,0)*r(2,0)*(0.25*a*r2 - 0.5) - 0.5*r(2,0)*r(2,0) + r(2,1)*r(2,1)*(0.25*a*r2 - 0.5) - 0.5*r(2,1)*r(2,1) + r(2,2)*r(2,2)*(0.25*a*r2 - 0.5) - 0.5*r(2,2)*r(2,2))*exp(a*(8.5*r2 + 8.0*r3)))*exp(-9.0*a*(r2 + r3))/(r2*r2*r2*r3*r3*r3);


	mat dujmat  = zeros<mat>( 4,3);
	dujmat(0,0) = a*r(0,0)*((-0.25*a*r0 + 1.0)*exp(a*(r0 + 0.5*r1)) + (0.5*a*r1 - 1.0)*exp(a*(0.5*r0 + r1)))*exp(-1.5*a*(r0 + r1))/r0;

	dujmat(0,1) = a*r(0,1)*((-0.25*a*r0 + 1.0)*exp(a*(r0 + 0.5*r1)) + (0.5*a*r1 - 1.0)*exp(a*(0.5*r0 + r1)))*exp(-1.5*a*(r0 + r1))/r0;

	dujmat(0,2) = a*r(0,2)*((-0.25*a*r0 + 1.0)*exp(a*(r0 + 0.5*r1)) + (0.5*a*r1 - 1.0)*exp(a*(0.5*r0 + r1)))*exp(-1.5*a*(r0 + r1))/r0;

	dujmat(1,0) = a*r(1,0)*(-(0.5*a*r0 - 1.0)*exp(a*(r0 + 0.5*r1)) + (0.25*a*r1 - 1.0)*exp(a*(0.5*r0 + r1)))*exp(-1.5*a*(r0 + r1))/r1;

	dujmat(1,1) = a*r(1,1)*(-(0.5*a*r0 - 1.0)*exp(a*(r0 + 0.5*r1)) + (0.25*a*r1 - 1.0)*exp(a*(0.5*r0 + r1)))*exp(-1.5*a*(r0 + r1))/r1;

	dujmat(1,2) = a*r(1,2)*(-(0.5*a*r0 - 1.0)*exp(a*(r0 + 0.5*r1)) + (0.25*a*r1 - 1.0)*exp(a*(0.5*r0 + r1)))*exp(-1.5*a*(r0 + r1))/r1;

	dujmat(2,0) = 0;

	dujmat(2,1) = 0;

	dujmat(2,2) = 0;

	dujmat(3,0) = 0;

	dujmat(3,1) = 0;

	dujmat(3,2) = 0;


	mat ddjmat  = zeros<mat>( 4,3);
	ddjmat(0,0) = 0;

	ddjmat(0,1) = 0;

	ddjmat(0,2) = 0;

	ddjmat(1,0) = 0;

	ddjmat(1,1) = 0;

	ddjmat(1,2) = 0;

	ddjmat(2,0) = a*r(2,0)*((-0.25*a*r2 + 1.0)*exp(a*(r2 + 0.5*r3)) + (0.5*a*r3 - 1.0)*exp(a*(0.5*r2 + r3)))*exp(-1.5*a*(r2 + r3))/r2;

	ddjmat(2,1) = a*r(2,1)*((-0.25*a*r2 + 1.0)*exp(a*(r2 + 0.5*r3)) + (0.5*a*r3 - 1.0)*exp(a*(0.5*r2 + r3)))*exp(-1.5*a*(r2 + r3))/r2;

	ddjmat(2,2) = a*r(2,2)*((-0.25*a*r2 + 1.0)*exp(a*(r2 + 0.5*r3)) + (0.5*a*r3 - 1.0)*exp(a*(0.5*r2 + r3)))*exp(-1.5*a*(r2 + r3))/r2;

	ddjmat(3,0) = a*r(3,0)*(-(0.5*a*r2 - 1.0)*exp(a*(r2 + 0.5*r3)) + (0.25*a*r3 - 1.0)*exp(a*(0.5*r2 + r3)))*exp(-1.5*a*(r2 + r3))/r3;

	ddjmat(3,1) = a*r(3,1)*(-(0.5*a*r2 - 1.0)*exp(a*(r2 + 0.5*r3)) + (0.25*a*r3 - 1.0)*exp(a*(0.5*r2 + r3)))*exp(-1.5*a*(r2 + r3))/r3;

	ddjmat(3,2) = a*r(3,2)*(-(0.5*a*r2 - 1.0)*exp(a*(r2 + 0.5*r3)) + (0.25*a*r3 - 1.0)*exp(a*(0.5*r2 + r3)))*exp(-1.5*a*(r2 + r3))/r3;

	double updet = -(-0.5*a*r0 + 1.0)*exp(-0.5*a*r0)*exp(-a*r1) + (-0.5*a*r1 + 1.0)*exp(-a*r0)*exp(-0.5*a*r1);

	double downdet = -(-0.5*a*r2 + 1.0)*exp(-0.5*a*r2)*exp(-a*r3) + (-0.5*a*r3 + 1.0)*exp(-a*r2)*exp(-0.5*a*r3);

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

	double dujgradient = 0.0;
	double ddjgradient = 0.0;

	for(int i = 0 ; i < nParticles; i++) {
		dujgradient += dot(dujmat.row(i),gradient.row(i));
		ddjgradient += dot(ddjmat.row(i),gradient.row(i));
	}

	return -0.5*d2u/updet - 0.5*d2d/downdet - 0.5*laplacian - dujgradient/updet - ddjgradient/downdet + coreEl + electronEl;
}
