#include <armadillo>
#include <math.h>
#include "Hydrogenic.h"
using namespace arma;
using namespace std;


Hydrogenic::Hydrogenic(int nDim,int nPart,int ch,double a,double b) {
	nDimensions = nDim;
	nParticles = nPart;
	charge = ch;
	alpha = a;
	beta = b;

}

void Hydrogenic::setR(double R){
 
}

double Hydrogenic::waveFunction(const mat &r,int nParticle,int orbital) {

	int i = nParticle;
	double rSingleParticle = 0;
	double argument = 0.0;

	for(int j = 0; j < nDimensions; j++) {
		rSingleParticle += r(i,j) * r(i,j);
	}

	argument += sqrt(rSingleParticle);
	switch (orbital) {
	case 0:
		return exp(-alpha*argument);
	case 1:
		return (-alpha*argument + 2)*exp(-1.0L/2.0L*alpha*argument);
	case 2:
		return r(i,1)*exp(-1.0L/2.0L*alpha*argument);
	case 3:
		return r(i,2)*exp(-1.0L/2.0L*alpha*argument);
	case 4:
		return r(i,0)*exp(-1.0L/2.0L*alpha*argument);
	}
}

rowvec Hydrogenic::gradient(const mat &r,int nParticle,int orbital) {

	int i = nParticle;
	double rSingleParticle = 0;
	double argument = 0.0;

	for(int j = 0; j < nDimensions; j++) {
		rSingleParticle += r(i,j) * r(i,j);
	}

	argument += sqrt(rSingleParticle);
	rowvec retvec(3);

	switch (orbital) {
	case 0:
		retvec(0) = -alpha*r(i,0)*exp(-alpha*argument)/argument;
		retvec(1) = -alpha*r(i,1)*exp(-alpha*argument)/argument;
		retvec(2) = -alpha*r(i,2)*exp(-alpha*argument)/argument;
		return retvec*(1);

	case 1:
		retvec(0) = (1.0L/2.0L)*alpha*r(i,0)*(alpha*argument - 4)*exp(-1.0L/2.0L*alpha*argument)/argument;
		retvec(1) = (1.0L/2.0L)*alpha*r(i,1)*(alpha*argument - 4)*exp(-1.0L/2.0L*alpha*argument)/argument;
		retvec(2) = (1.0L/2.0L)*alpha*r(i,2)*(alpha*argument - 4)*exp(-1.0L/2.0L*alpha*argument)/argument;
		return retvec*(1);

	case 2:
		retvec(0) = -1.0L/2.0L*alpha*r(i,0)*r(i,1)*exp(-1.0L/2.0L*alpha*argument)/argument;
		retvec(1) = (-1.0L/2.0L*alpha*(r(i,1)*r(i,1)) + argument)*exp(-1.0L/2.0L*alpha*argument)/argument;
		retvec(2) = -1.0L/2.0L*alpha*r(i,1)*r(i,2)*exp(-1.0L/2.0L*alpha*argument)/argument;
		return retvec*(1);

	case 3:
		retvec(0) = -1.0L/2.0L*alpha*r(i,0)*r(i,2)*exp(-1.0L/2.0L*alpha*argument)/argument;
		retvec(1) = -1.0L/2.0L*alpha*r(i,1)*r(i,2)*exp(-1.0L/2.0L*alpha*argument)/argument;
		retvec(2) = (-1.0L/2.0L*alpha*(r(i,2)*r(i,2)) + argument)*exp(-1.0L/2.0L*alpha*argument)/argument;
		return retvec*(1);

	case 4:
		retvec(0) = (-1.0L/2.0L*alpha*(r(i,0)*r(i,0)) + argument)*exp(-1.0L/2.0L*alpha*argument)/argument;
		retvec(1) = -1.0L/2.0L*alpha*r(i,0)*r(i,1)*exp(-1.0L/2.0L*alpha*argument)/argument;
		retvec(2) = -1.0L/2.0L*alpha*r(i,0)*r(i,2)*exp(-1.0L/2.0L*alpha*argument)/argument;
		return retvec*(1);

	}
}
double Hydrogenic::laplacian(const mat &r,int nParticle,int orbital) {

	int i = nParticle;
	double rSingleParticle = 0;
	double argument = 0.0;

	for(int j = 0; j < nDimensions; j++) {
		rSingleParticle += r(i,j) * r(i,j);
	}

	argument += sqrt(rSingleParticle);
	rowvec retvec(3);

	switch (orbital) {
	case 0:
		return ((alpha*alpha)*(r(i,0)*r(i,0))/(argument*argument) + (alpha*alpha)*(r(i,1)*r(i,1))/(argument*argument) + (alpha*alpha)*(r(i,2)*r(i,2))/(argument*argument) - 3*alpha/argument + alpha*(r(i,0)*r(i,0))/(argument*argument*argument) + alpha*(r(i,1)*r(i,1))/(argument*argument*argument) + alpha*(r(i,2)*r(i,2))/(argument*argument*argument))*exp(-alpha*argument);
	case 1:
		return ((1.0L/4.0L)*(alpha*alpha)*(r(i,0)*r(i,0))*(-alpha*argument + 2)/(argument*argument) + (alpha*alpha)*(r(i,0)*r(i,0))/(argument*argument) + (1.0L/4.0L)*(alpha*alpha)*(r(i,1)*r(i,1))*(-alpha*argument + 2)/(argument*argument) + (alpha*alpha)*(r(i,1)*r(i,1))/(argument*argument) + (1.0L/4.0L)*(alpha*alpha)*(r(i,2)*r(i,2))*(-alpha*argument + 2)/(argument*argument) + (alpha*alpha)*(r(i,2)*r(i,2))/(argument*argument) - 3.0L/2.0L*alpha*(-alpha*argument + 2)/argument - 3*alpha/argument + (1.0L/2.0L)*alpha*(r(i,0)*r(i,0))*(-alpha*argument + 2)/(argument*argument*argument) + alpha*(r(i,0)*r(i,0))/(argument*argument*argument) + (1.0L/2.0L)*alpha*(r(i,1)*r(i,1))*(-alpha*argument + 2)/(argument*argument*argument) + alpha*(r(i,1)*r(i,1))/(argument*argument*argument) + (1.0L/2.0L)*alpha*(r(i,2)*r(i,2))*(-alpha*argument + 2)/(argument*argument*argument) + alpha*(r(i,2)*r(i,2))/(argument*argument*argument))/sqrt(exp(alpha*argument));
	case 2:
		return ((1.0L/4.0L)*(alpha*alpha)*(r(i,0)*r(i,0))*r(i,1)/(argument*argument) + (1.0L/4.0L)*(alpha*alpha)*(r(i,1)*r(i,1)*r(i,1))/(argument*argument) + (1.0L/4.0L)*(alpha*alpha)*r(i,1)*(r(i,2)*r(i,2))/(argument*argument) - 5.0L/2.0L*alpha*r(i,1)/argument + (1.0L/2.0L)*alpha*(r(i,0)*r(i,0))*r(i,1)/(argument*argument*argument) + (1.0L/2.0L)*alpha*(r(i,1)*r(i,1)*r(i,1))/(argument*argument*argument) + (1.0L/2.0L)*alpha*r(i,1)*(r(i,2)*r(i,2))/(argument*argument*argument))/sqrt(exp(alpha*argument));
	case 3:
		return ((1.0L/4.0L)*(alpha*alpha)*(r(i,0)*r(i,0))*r(i,2)/(argument*argument) + (1.0L/4.0L)*(alpha*alpha)*(r(i,1)*r(i,1))*r(i,2)/(argument*argument) + (1.0L/4.0L)*(alpha*alpha)*(r(i,2)*r(i,2)*r(i,2))/(argument*argument) - 5.0L/2.0L*alpha*r(i,2)/argument + (1.0L/2.0L)*alpha*(r(i,0)*r(i,0))*r(i,2)/(argument*argument*argument) + (1.0L/2.0L)*alpha*(r(i,1)*r(i,1))*r(i,2)/(argument*argument*argument) + (1.0L/2.0L)*alpha*(r(i,2)*r(i,2)*r(i,2))/(argument*argument*argument))/sqrt(exp(alpha*argument));
	case 4:
		return ((1.0L/4.0L)*(alpha*alpha)*(r(i,0)*r(i,0)*r(i,0))/(argument*argument) + (1.0L/4.0L)*(alpha*alpha)*r(i,0)*(r(i,1)*r(i,1))/(argument*argument) + (1.0L/4.0L)*(alpha*alpha)*r(i,0)*(r(i,2)*r(i,2))/(argument*argument) - 5.0L/2.0L*alpha*r(i,0)/argument + (1.0L/2.0L)*alpha*(r(i,0)*r(i,0)*r(i,0))/(argument*argument*argument) + (1.0L/2.0L)*alpha*r(i,0)*(r(i,1)*r(i,1))/(argument*argument*argument) + (1.0L/2.0L)*alpha*r(i,0)*(r(i,2)*r(i,2))/(argument*argument*argument))/sqrt(exp(alpha*argument));
	}

}
