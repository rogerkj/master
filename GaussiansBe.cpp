#include <armadillo>
#include <math.h>
#include "Gaussians.h"
using namespace arma;
using namespace std;


void Gaussians::setR(double R){
}

Gaussians::Gaussians(int nDim,int nPart,int ch,double a,double b) {
	nDimensions = nDim;
	nParticles = nPart;
	charge = ch;
	alpha = a;
	beta = b;

}


double Gaussians::waveFunction(const mat &r,int nParticle,int orbital) {

	int i = nParticle;
	double rSingleParticle = 0;
	double argument = 0.0;

	for(int j = 0; j < nDimensions; j++) {
		rSingleParticle += r(i,j) * r(i,j);
	}

	switch (orbital) {
	case 0:
		return 0.440630314098264*exp(-6.36242139*rSingleParticle) + 0.426158447617189*exp(-1.158923*rSingleParticle) + 0.132814973357573*exp(-0.31364979*rSingleParticle);
	}
}

rowvec Gaussians::gradient(const mat &r,int nParticle,int orbital) {

	int i = nParticle;
	double rSingleParticle = 0;
	double argument = 0.0;

	for(int j = 0; j < nDimensions; j++) {
		rSingleParticle += r(i,j) * r(i,j);
	}

	rowvec retvec(3);

	switch (orbital) {
	case 0:
		retvec(0) = -5.60695147100242*r(i,0)*exp(-6.36242139*rSingleParticle) - 0.987769653175711*r(i,0)*exp(-1.158923*rSingleParticle) - 0.083314777004917*r(i,0)*exp(-0.31364979*rSingleParticle);
		retvec(1) = -5.60695147100242*r(i,1)*exp(-6.36242139*rSingleParticle) - 0.987769653175711*r(i,1)*exp(-1.158923*rSingleParticle) - 0.083314777004917*r(i,1)*exp(-0.31364979*rSingleParticle);
		retvec(2) = -5.60695147100242*r(i,2)*exp(-6.36242139*rSingleParticle) - 0.987769653175711*r(i,2)*exp(-1.158923*rSingleParticle) - 0.083314777004917*r(i,2)*exp(-0.31364979*rSingleParticle);
		return retvec;

	}
}
double Gaussians::laplacian(const mat &r,int nParticle,int orbital) {

	int i = nParticle;
	double rSingleParticle = 0;
	double argument = 0.0;

	for(int j = 0; j < nDimensions; j++) {
		rSingleParticle += r(i,j) * r(i,j);
	}

	rowvec retvec(3);

	switch (orbital) {
	case 0:
		return (0.0522633246229781*(r(i,0)*r(i,0)) + 0.0522633246229781*(r(i,1)*r(i,1)) + 0.0522633246229781*(r(i,2)*r(i,2)) - 0.249944331014751)*exp(-0.31364979*rSingleParticle) + (2.28949793953471*(r(i,0)*r(i,0)) + 2.28949793953471*(r(i,1)*r(i,1)) + 2.28949793953471*(r(i,2)*r(i,2)) - 2.96330895952713)*exp(-1.158923*rSingleParticle) + (71.3475759435955*(r(i,0)*r(i,0)) + 71.3475759435955*(r(i,1)*r(i,1)) + 71.3475759435955*(r(i,2)*r(i,2)) - 16.8208544130073)*exp(-6.36242139*rSingleParticle);
	}

}
