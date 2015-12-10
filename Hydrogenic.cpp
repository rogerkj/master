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


void Hydrogenic::setR(double R){} 

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
	case 5:
		return (2*(alpha*alpha)*((r(i,0)*r(i,0)) + (r(i,1)*r(i,1)) + (r(i,2)*r(i,2))) - 18*alpha*argument + 27)*exp(-1.0L/3.0L*alpha*argument);
	case 6:
		return r(i,1)*(alpha*argument - 6)*exp(-1.0L/3.0L*alpha*argument);
	case 7:
		return r(i,2)*(alpha*argument - 6)*exp(-1.0L/3.0L*alpha*argument);
	case 8:
		return r(i,0)*(alpha*argument - 6)*exp(-1.0L/3.0L*alpha*argument);
	case 9:
		return r(i,0)*r(i,1)*exp(-1.0L/3.0L*alpha*argument);
	case 10:
		return r(i,1)*r(i,2)*exp(-1.0L/3.0L*alpha*argument);
	case 11:
		return ((r(i,0)*r(i,0)) + (r(i,1)*r(i,1)) - 2*(r(i,2)*r(i,2)))*exp(-1.0L/3.0L*alpha*argument);
	case 12:
		return r(i,0)*r(i,2)*exp(-1.0L/3.0L*alpha*argument);
	case 13:
		return ((r(i,0)*r(i,0)) - (r(i,1)*r(i,1)))*exp(-1.0L/3.0L*alpha*argument);
	case 14:
		return ((alpha*alpha)*(r(i,0)*r(i,0))*(-alpha*argument + 24) + (alpha*alpha)*(r(i,1)*r(i,1))*(-alpha*argument + 24) + (alpha*alpha)*(r(i,2)*r(i,2))*(-alpha*argument + 24) - 144*alpha*argument + 192)*exp(-1.0L/4.0L*alpha*argument);
	case 15:
		return r(i,1)*((alpha*alpha)*rSingleParticle - 20*alpha*argument + 80)*exp(-1.0L/4.0L*alpha*argument);
	case 16:
		return r(i,2)*((alpha*alpha)*rSingleParticle - 20*alpha*argument + 80)*exp(-1.0L/4.0L*alpha*argument);
	case 17:
		return r(i,0)*((alpha*alpha)*rSingleParticle - 20*alpha*argument + 80)*exp(-1.0L/4.0L*alpha*argument);
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
		retvec(0) = r(i,0);
		retvec(1) = r(i,1);
		retvec(2) = r(i,2);
		return retvec*(-alpha*exp(-alpha*argument)/argument);

	
	
	case 1:
		retvec(0) = r(i,0);
		retvec(1) = r(i,1);
		retvec(2) = r(i,2);
		return retvec*((1.0L/2.0L)*alpha*(alpha*argument - 4)*exp(-1.0L/2.0L*alpha*argument)/argument);

	
	
	case 2:
		retvec(0) = -1.0L/2.0L*alpha*r(i,0)*r(i,1);
		retvec(1) = -1.0L/2.0L*alpha*(r(i,1)*r(i,1)) + argument;
		retvec(2) = -1.0L/2.0L*alpha*r(i,1)*r(i,2);
		return retvec*(exp(-1.0L/2.0L*alpha*argument)/argument);

	
	
	case 3:
		retvec(0) = -1.0L/2.0L*alpha*r(i,0)*r(i,2);
		retvec(1) = -1.0L/2.0L*alpha*r(i,1)*r(i,2);
		retvec(2) = -1.0L/2.0L*alpha*(r(i,2)*r(i,2)) + argument;
		return retvec*(exp(-1.0L/2.0L*alpha*argument)/argument);

	
	
	case 4:
		retvec(0) = -1.0L/2.0L*alpha*(r(i,0)*r(i,0)) + argument;
		retvec(1) = -1.0L/2.0L*alpha*r(i,0)*r(i,1);
		retvec(2) = -1.0L/2.0L*alpha*r(i,0)*r(i,2);
		return retvec*(exp(-1.0L/2.0L*alpha*argument)/argument);

	
	
	case 5:
		retvec(0) = r(i,0);
		retvec(1) = r(i,1);
		retvec(2) = r(i,2);
		return retvec*(-1.0L/3.0L*alpha*(2*(alpha*alpha)*((r(i,0)*r(i,0)) + (r(i,1)*r(i,1)) + (r(i,2)*r(i,2))) - 30*alpha*argument + 81)*exp(-1.0L/3.0L*alpha*argument)/argument);

	
	
	case 6:
		retvec(0) = (1.0L/3.0L)*alpha*r(i,0)*r(i,1)*(-alpha*argument + 9);
		retvec(1) = -1.0L/3.0L*alpha*(r(i,1)*r(i,1))*(alpha*argument - 9) + argument*(alpha*argument - 6);
		retvec(2) = (1.0L/3.0L)*alpha*r(i,1)*r(i,2)*(-alpha*argument + 9);
		return retvec*(exp(-1.0L/3.0L*alpha*argument)/argument);

	
	
	case 7:
		retvec(0) = (1.0L/3.0L)*alpha*r(i,0)*r(i,2)*(-alpha*argument + 9);
		retvec(1) = (1.0L/3.0L)*alpha*r(i,1)*r(i,2)*(-alpha*argument + 9);
		retvec(2) = -1.0L/3.0L*alpha*(r(i,2)*r(i,2))*(alpha*argument - 9) + argument*(alpha*argument - 6);
		return retvec*(exp(-1.0L/3.0L*alpha*argument)/argument);

	
	
	case 8:
		retvec(0) = -1.0L/3.0L*alpha*(r(i,0)*r(i,0))*(alpha*argument - 9) + argument*(alpha*argument - 6);
		retvec(1) = (1.0L/3.0L)*alpha*r(i,0)*r(i,1)*(-alpha*argument + 9);
		retvec(2) = (1.0L/3.0L)*alpha*r(i,0)*r(i,2)*(-alpha*argument + 9);
		return retvec*(exp(-1.0L/3.0L*alpha*argument)/argument);

	
	
	case 9:
		retvec(0) = (1.0L/3.0L)*r(i,1)*(-alpha*(r(i,0)*r(i,0)) + 3*argument);
		retvec(1) = (1.0L/3.0L)*r(i,0)*(-alpha*(r(i,1)*r(i,1)) + 3*argument);
		retvec(2) = -1.0L/3.0L*alpha*r(i,0)*r(i,1)*r(i,2);
		return retvec*(exp(-1.0L/3.0L*alpha*argument)/argument);

	
	
	case 10:
		retvec(0) = -1.0L/3.0L*alpha*r(i,0)*r(i,1)*r(i,2);
		retvec(1) = (1.0L/3.0L)*r(i,2)*(-alpha*(r(i,1)*r(i,1)) + 3*argument);
		retvec(2) = (1.0L/3.0L)*r(i,1)*(-alpha*(r(i,2)*r(i,2)) + 3*argument);
		return retvec*(exp(-1.0L/3.0L*alpha*argument)/argument);

	
	
	case 11:
		retvec(0) = -1.0L/3.0L*r(i,0)*(alpha*((r(i,0)*r(i,0)) + (r(i,1)*r(i,1)) - 2*(r(i,2)*r(i,2))) - 6*argument);
		retvec(1) = -1.0L/3.0L*r(i,1)*(alpha*((r(i,0)*r(i,0)) + (r(i,1)*r(i,1)) - 2*(r(i,2)*r(i,2))) - 6*argument);
		retvec(2) = -1.0L/3.0L*r(i,2)*(alpha*((r(i,0)*r(i,0)) + (r(i,1)*r(i,1)) - 2*(r(i,2)*r(i,2))) + 12*argument);
		return retvec*(exp(-1.0L/3.0L*alpha*argument)/argument);

	
	
	case 12:
		retvec(0) = (1.0L/3.0L)*r(i,2)*(-alpha*(r(i,0)*r(i,0)) + 3*argument);
		retvec(1) = -1.0L/3.0L*alpha*r(i,0)*r(i,1)*r(i,2);
		retvec(2) = (1.0L/3.0L)*r(i,0)*(-alpha*(r(i,2)*r(i,2)) + 3*argument);
		return retvec*(exp(-1.0L/3.0L*alpha*argument)/argument);

	
	
	case 13:
		retvec(0) = -1.0L/3.0L*r(i,0)*(alpha*((r(i,0)*r(i,0)) - (r(i,1)*r(i,1))) - 6*argument);
		retvec(1) = -1.0L/3.0L*r(i,1)*(alpha*((r(i,0)*r(i,0)) - (r(i,1)*r(i,1))) + 6*argument);
		retvec(2) = (1.0L/3.0L)*alpha*r(i,2)*(-(r(i,0)*r(i,0)) + (r(i,1)*r(i,1)));
		return retvec*(exp(-1.0L/3.0L*alpha*argument)/argument);

	
	
	case 14:
		retvec(0) = r(i,0);
		retvec(1) = r(i,1);
		retvec(2) = r(i,2);
		return retvec*((1.0L/4.0L)*alpha*((alpha*alpha*alpha)*argument*(r(i,0)*r(i,0)) + (alpha*alpha*alpha)*argument*(r(i,1)*r(i,1)) + (alpha*alpha*alpha)*argument*(r(i,2)*r(i,2)) - 8*(alpha*alpha)*(argument*argument) - 28*(alpha*alpha)*(r(i,0)*r(i,0)) - 28*(alpha*alpha)*(r(i,1)*r(i,1)) - 28*(alpha*alpha)*(r(i,2)*r(i,2)) + 336*alpha*argument - 768)*exp(-1.0L/4.0L*alpha*argument)/argument);

	
	
	case 15:
		retvec(0) = (1.0L/4.0L)*alpha*r(i,0)*r(i,1)*(-(alpha*alpha)*rSingleParticle + 28*alpha*argument - 160);
		retvec(1) = -1.0L/4.0L*alpha*(r(i,1)*r(i,1))*((alpha*alpha)*rSingleParticle - 28*alpha*argument + 160) + argument*((alpha*alpha)*rSingleParticle - 20*alpha*argument + 80);
		retvec(2) = (1.0L/4.0L)*alpha*r(i,1)*r(i,2)*(-(alpha*alpha)*rSingleParticle + 28*alpha*argument - 160);
		return retvec*(exp(-1.0L/4.0L*alpha*argument)/argument);

	
	
	case 16:
		retvec(0) = (1.0L/4.0L)*alpha*r(i,0)*r(i,2)*(-(alpha*alpha)*rSingleParticle + 28*alpha*argument - 160);
		retvec(1) = (1.0L/4.0L)*alpha*r(i,1)*r(i,2)*(-(alpha*alpha)*rSingleParticle + 28*alpha*argument - 160);
		retvec(2) = -1.0L/4.0L*alpha*(r(i,2)*r(i,2))*((alpha*alpha)*rSingleParticle - 28*alpha*argument + 160) + argument*((alpha*alpha)*rSingleParticle - 20*alpha*argument + 80);
		return retvec*(exp(-1.0L/4.0L*alpha*argument)/argument);

	
	
	case 17:
		retvec(0) = -1.0L/4.0L*alpha*(r(i,0)*r(i,0))*((alpha*alpha)*rSingleParticle - 28*alpha*argument + 160) + argument*((alpha*alpha)*rSingleParticle - 20*alpha*argument + 80);
		retvec(1) = (1.0L/4.0L)*alpha*r(i,0)*r(i,1)*(-(alpha*alpha)*rSingleParticle + 28*alpha*argument - 160);
		retvec(2) = (1.0L/4.0L)*alpha*r(i,0)*r(i,2)*(-(alpha*alpha)*rSingleParticle + 28*alpha*argument - 160);
		return retvec*(exp(-1.0L/4.0L*alpha*argument)/argument);

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
		return alpha*(alpha*argument*((r(i,0)*r(i,0)) + (r(i,1)*r(i,1)) + (r(i,2)*r(i,2))) - 3*(argument*argument) + (r(i,0)*r(i,0)) + (r(i,1)*r(i,1)) + (r(i,2)*r(i,2)))*exp(-alpha*argument)/(argument*argument*argument);
	case 1:
		return (1.0L/4.0L)*alpha*(-(alpha*alpha)*(argument*argument)*(r(i,0)*r(i,0)) - (alpha*alpha)*(argument*argument)*(r(i,1)*r(i,1)) - (alpha*alpha)*(argument*argument)*(r(i,2)*r(i,2)) + 6*alpha*(argument*argument*argument) + 4*alpha*argument*(r(i,0)*r(i,0)) + 4*alpha*argument*(r(i,1)*r(i,1)) + 4*alpha*argument*(r(i,2)*r(i,2)) - 24*(argument*argument) + 8*(r(i,0)*r(i,0)) + 8*(r(i,1)*r(i,1)) + 8*(r(i,2)*r(i,2)))*exp(-1.0L/2.0L*alpha*argument)/(argument*argument*argument);
	case 2:
		return (1.0L/4.0L)*alpha*r(i,1)*(alpha*argument*((r(i,0)*r(i,0)) + (r(i,1)*r(i,1)) + (r(i,2)*r(i,2))) - 10*(argument*argument) + 2*(r(i,0)*r(i,0)) + 2*(r(i,1)*r(i,1)) + 2*(r(i,2)*r(i,2)))*exp(-1.0L/2.0L*alpha*argument)/(argument*argument*argument);
	case 3:
		return (1.0L/4.0L)*alpha*r(i,2)*(alpha*argument*((r(i,0)*r(i,0)) + (r(i,1)*r(i,1)) + (r(i,2)*r(i,2))) - 10*(argument*argument) + 2*(r(i,0)*r(i,0)) + 2*(r(i,1)*r(i,1)) + 2*(r(i,2)*r(i,2)))*exp(-1.0L/2.0L*alpha*argument)/(argument*argument*argument);
	case 4:
		return (1.0L/4.0L)*alpha*r(i,0)*(alpha*argument*((r(i,0)*r(i,0)) + (r(i,1)*r(i,1)) + (r(i,2)*r(i,2))) - 10*(argument*argument) + 2*(r(i,0)*r(i,0)) + 2*(r(i,1)*r(i,1)) + 2*(r(i,2)*r(i,2)))*exp(-1.0L/2.0L*alpha*argument)/(argument*argument*argument);
	case 5:
		return (1.0L/9.0L)*alpha*(108*alpha*(argument*argument*argument) + alpha*argument*(-12*(r(i,0)*r(i,0))*(2*alpha*argument - 9) + (r(i,0)*r(i,0))*(2*(alpha*alpha)*((r(i,0)*r(i,0)) + (r(i,1)*r(i,1)) + (r(i,2)*r(i,2))) - 18*alpha*argument + 27) - 12*(r(i,1)*r(i,1))*(2*alpha*argument - 9) + (r(i,1)*r(i,1))*(2*(alpha*alpha)*((r(i,0)*r(i,0)) + (r(i,1)*r(i,1)) + (r(i,2)*r(i,2))) - 18*alpha*argument + 27) - 12*(r(i,2)*r(i,2))*(2*alpha*argument - 9) + (r(i,2)*r(i,2))*(2*(alpha*alpha)*((r(i,0)*r(i,0)) + (r(i,1)*r(i,1)) + (r(i,2)*r(i,2))) - 18*alpha*argument + 27)) - 9*(argument*argument)*(2*(alpha*alpha)*((r(i,0)*r(i,0)) + (r(i,1)*r(i,1)) + (r(i,2)*r(i,2))) - 18*alpha*argument + 27) - 486*(argument*argument) + 162*(r(i,0)*r(i,0)) + 162*(r(i,1)*r(i,1)) + 162*(r(i,2)*r(i,2)) + 3*((r(i,0)*r(i,0)) + (r(i,1)*r(i,1)) + (r(i,2)*r(i,2)))*(2*(alpha*alpha)*((r(i,0)*r(i,0)) + (r(i,1)*r(i,1)) + (r(i,2)*r(i,2))) - 18*alpha*argument + 27))*exp(-1.0L/3.0L*alpha*argument)/(argument*argument*argument);
	case 6:
		return (1.0L/9.0L)*alpha*r(i,1)*((alpha*alpha)*(argument*argument)*(r(i,0)*r(i,0)) + (alpha*alpha)*(argument*argument)*(r(i,1)*r(i,1)) + (alpha*alpha)*(argument*argument)*(r(i,2)*r(i,2)) - 15*alpha*(argument*argument*argument) - 9*alpha*argument*(r(i,0)*r(i,0)) - 9*alpha*argument*(r(i,1)*r(i,1)) - 9*alpha*argument*(r(i,2)*r(i,2)) + 135*(argument*argument) - 27*(r(i,0)*r(i,0)) - 27*(r(i,1)*r(i,1)) - 27*(r(i,2)*r(i,2)))*exp(-1.0L/3.0L*alpha*argument)/(argument*argument*argument);
	case 7:
		return (1.0L/9.0L)*alpha*r(i,2)*((alpha*alpha)*(argument*argument)*(r(i,0)*r(i,0)) + (alpha*alpha)*(argument*argument)*(r(i,1)*r(i,1)) + (alpha*alpha)*(argument*argument)*(r(i,2)*r(i,2)) - 15*alpha*(argument*argument*argument) - 9*alpha*argument*(r(i,0)*r(i,0)) - 9*alpha*argument*(r(i,1)*r(i,1)) - 9*alpha*argument*(r(i,2)*r(i,2)) + 135*(argument*argument) - 27*(r(i,0)*r(i,0)) - 27*(r(i,1)*r(i,1)) - 27*(r(i,2)*r(i,2)))*exp(-1.0L/3.0L*alpha*argument)/(argument*argument*argument);
	case 8:
		return (1.0L/9.0L)*alpha*r(i,0)*((alpha*alpha)*(argument*argument)*(r(i,0)*r(i,0)) + (alpha*alpha)*(argument*argument)*(r(i,1)*r(i,1)) + (alpha*alpha)*(argument*argument)*(r(i,2)*r(i,2)) - 15*alpha*(argument*argument*argument) - 9*alpha*argument*(r(i,0)*r(i,0)) - 9*alpha*argument*(r(i,1)*r(i,1)) - 9*alpha*argument*(r(i,2)*r(i,2)) + 135*(argument*argument) - 27*(r(i,0)*r(i,0)) - 27*(r(i,1)*r(i,1)) - 27*(r(i,2)*r(i,2)))*exp(-1.0L/3.0L*alpha*argument)/(argument*argument*argument);
	case 9:
		return (1.0L/9.0L)*alpha*r(i,0)*r(i,1)*(alpha*argument*((r(i,0)*r(i,0)) + (r(i,1)*r(i,1)) + (r(i,2)*r(i,2))) - 21*(argument*argument) + 3*(r(i,0)*r(i,0)) + 3*(r(i,1)*r(i,1)) + 3*(r(i,2)*r(i,2)))*exp(-1.0L/3.0L*alpha*argument)/(argument*argument*argument);
	case 10:
		return (1.0L/9.0L)*alpha*r(i,1)*r(i,2)*(alpha*argument*((r(i,0)*r(i,0)) + (r(i,1)*r(i,1)) + (r(i,2)*r(i,2))) - 21*(argument*argument) + 3*(r(i,0)*r(i,0)) + 3*(r(i,1)*r(i,1)) + 3*(r(i,2)*r(i,2)))*exp(-1.0L/3.0L*alpha*argument)/(argument*argument*argument);
	case 11:
		return (1.0L/9.0L)*alpha*(alpha*argument*((r(i,0)*r(i,0)) + (r(i,1)*r(i,1)) - 2*(r(i,2)*r(i,2)))*((r(i,0)*r(i,0)) + (r(i,1)*r(i,1)) + (r(i,2)*r(i,2))) + 21*(argument*argument)*(-(r(i,0)*r(i,0)) - (r(i,1)*r(i,1)) + 2*(r(i,2)*r(i,2))) + 3*((r(i,0)*r(i,0)) + (r(i,1)*r(i,1)) - 2*(r(i,2)*r(i,2)))*((r(i,0)*r(i,0)) + (r(i,1)*r(i,1)) + (r(i,2)*r(i,2))))*exp(-1.0L/3.0L*alpha*argument)/(argument*argument*argument);
	case 12:
		return (1.0L/9.0L)*alpha*r(i,0)*r(i,2)*(alpha*argument*((r(i,0)*r(i,0)) + (r(i,1)*r(i,1)) + (r(i,2)*r(i,2))) - 21*(argument*argument) + 3*(r(i,0)*r(i,0)) + 3*(r(i,1)*r(i,1)) + 3*(r(i,2)*r(i,2)))*exp(-1.0L/3.0L*alpha*argument)/(argument*argument*argument);
	case 13:
		return (1.0L/9.0L)*alpha*(alpha*argument*((r(i,0)*r(i,0)) - (r(i,1)*r(i,1)))*((r(i,0)*r(i,0)) + (r(i,1)*r(i,1)) + (r(i,2)*r(i,2))) + 21*(argument*argument)*(-(r(i,0)*r(i,0)) + (r(i,1)*r(i,1))) + 3*((r(i,0)*r(i,0)) - (r(i,1)*r(i,1)))*((r(i,0)*r(i,0)) + (r(i,1)*r(i,1)) + (r(i,2)*r(i,2))))*exp(-1.0L/3.0L*alpha*argument)/(argument*argument*argument);
	case 14:
		return (1.0L/16.0L)*alpha*(-(alpha*alpha*alpha*alpha)*(argument*argument)*(r(i,0)*r(i,0)*r(i,0)*r(i,0)) - 2*(alpha*alpha*alpha*alpha)*(argument*argument)*(r(i,0)*r(i,0))*(r(i,1)*r(i,1)) - 2*(alpha*alpha*alpha*alpha)*(argument*argument)*(r(i,0)*r(i,0))*(r(i,2)*r(i,2)) - (alpha*alpha*alpha*alpha)*(argument*argument)*(r(i,1)*r(i,1)*r(i,1)*r(i,1)) - 2*(alpha*alpha*alpha*alpha)*(argument*argument)*(r(i,1)*r(i,1))*(r(i,2)*r(i,2)) - (alpha*alpha*alpha*alpha)*(argument*argument)*(r(i,2)*r(i,2)*r(i,2)*r(i,2)) + 28*(alpha*alpha*alpha)*(argument*argument*argument)*(r(i,0)*r(i,0)) + 28*(alpha*alpha*alpha)*(argument*argument*argument)*(r(i,1)*r(i,1)) + 28*(alpha*alpha*alpha)*(argument*argument*argument)*(r(i,2)*r(i,2)) + 28*(alpha*alpha*alpha)*argument*(r(i,0)*r(i,0)*r(i,0)*r(i,0)) + 56*(alpha*alpha*alpha)*argument*(r(i,0)*r(i,0))*(r(i,1)*r(i,1)) + 56*(alpha*alpha*alpha)*argument*(r(i,0)*r(i,0))*(r(i,2)*r(i,2)) + 28*(alpha*alpha*alpha)*argument*(r(i,1)*r(i,1)*r(i,1)*r(i,1)) + 56*(alpha*alpha*alpha)*argument*(r(i,1)*r(i,1))*(r(i,2)*r(i,2)) + 28*(alpha*alpha*alpha)*argument*(r(i,2)*r(i,2)*r(i,2)*r(i,2)) - 96*(alpha*alpha)*(argument*argument*argument*argument) - 928*(alpha*alpha)*(argument*argument)*(r(i,0)*r(i,0)) - 928*(alpha*alpha)*(argument*argument)*(r(i,1)*r(i,1)) - 928*(alpha*alpha)*(argument*argument)*(r(i,2)*r(i,2)) + 112*(alpha*alpha)*(r(i,0)*r(i,0)*r(i,0)*r(i,0)) + 224*(alpha*alpha)*(r(i,0)*r(i,0))*(r(i,1)*r(i,1)) + 224*(alpha*alpha)*(r(i,0)*r(i,0))*(r(i,2)*r(i,2)) + 112*(alpha*alpha)*(r(i,1)*r(i,1)*r(i,1)*r(i,1)) + 224*(alpha*alpha)*(r(i,1)*r(i,1))*(r(i,2)*r(i,2)) + 112*(alpha*alpha)*(r(i,2)*r(i,2)*r(i,2)*r(i,2)) + 4032*alpha*(argument*argument*argument) + 768*alpha*argument*(r(i,0)*r(i,0)) + 768*alpha*argument*(r(i,1)*r(i,1)) + 768*alpha*argument*(r(i,2)*r(i,2)) - 9216*(argument*argument) + 3072*(r(i,0)*r(i,0)) + 3072*(r(i,1)*r(i,1)) + 3072*(r(i,2)*r(i,2)))*exp(-1.0L/4.0L*alpha*argument)/(argument*argument*argument);
	case 15:
		return (1.0L/16.0L)*alpha*r(i,1)*(96*alpha*(argument*argument*argument) + alpha*argument*((r(i,0)*r(i,0))*(-16*alpha*argument + 160) + (r(i,0)*r(i,0))*((alpha*alpha)*rSingleParticle - 20*alpha*argument + 80) + (r(i,1)*r(i,1))*(-16*alpha*argument + 160) + (r(i,1)*r(i,1))*((alpha*alpha)*rSingleParticle - 20*alpha*argument + 80) + (r(i,2)*r(i,2))*(-16*alpha*argument + 160) + (r(i,2)*r(i,2))*((alpha*alpha)*rSingleParticle - 20*alpha*argument + 80)) + (argument*argument)*(-20*(alpha*alpha)*rSingleParticle + 464*alpha*argument - 2240) - 960*(argument*argument) + 320*(r(i,0)*r(i,0)) + 320*(r(i,1)*r(i,1)) + 320*(r(i,2)*r(i,2)) + 4*((r(i,0)*r(i,0)) + (r(i,1)*r(i,1)) + (r(i,2)*r(i,2)))*((alpha*alpha)*rSingleParticle - 20*alpha*argument + 80))*exp(-1.0L/4.0L*alpha*argument)/(argument*argument*argument);
	case 16:
		return (1.0L/16.0L)*alpha*r(i,2)*(96*alpha*(argument*argument*argument) + alpha*argument*((r(i,0)*r(i,0))*(-16*alpha*argument + 160) + (r(i,0)*r(i,0))*((alpha*alpha)*rSingleParticle - 20*alpha*argument + 80) + (r(i,1)*r(i,1))*(-16*alpha*argument + 160) + (r(i,1)*r(i,1))*((alpha*alpha)*rSingleParticle - 20*alpha*argument + 80) + (r(i,2)*r(i,2))*(-16*alpha*argument + 160) + (r(i,2)*r(i,2))*((alpha*alpha)*rSingleParticle - 20*alpha*argument + 80)) + (argument*argument)*(-20*(alpha*alpha)*rSingleParticle + 464*alpha*argument - 2240) - 960*(argument*argument) + 320*(r(i,0)*r(i,0)) + 320*(r(i,1)*r(i,1)) + 320*(r(i,2)*r(i,2)) + 4*((r(i,0)*r(i,0)) + (r(i,1)*r(i,1)) + (r(i,2)*r(i,2)))*((alpha*alpha)*rSingleParticle - 20*alpha*argument + 80))*exp(-1.0L/4.0L*alpha*argument)/(argument*argument*argument);
	case 17:
		return (1.0L/16.0L)*alpha*r(i,0)*(96*alpha*(argument*argument*argument) + alpha*argument*((r(i,0)*r(i,0))*(-16*alpha*argument + 160) + (r(i,0)*r(i,0))*((alpha*alpha)*rSingleParticle - 20*alpha*argument + 80) + (r(i,1)*r(i,1))*(-16*alpha*argument + 160) + (r(i,1)*r(i,1))*((alpha*alpha)*rSingleParticle - 20*alpha*argument + 80) + (r(i,2)*r(i,2))*(-16*alpha*argument + 160) + (r(i,2)*r(i,2))*((alpha*alpha)*rSingleParticle - 20*alpha*argument + 80)) + (argument*argument)*(-20*(alpha*alpha)*rSingleParticle + 464*alpha*argument - 2240) - 960*(argument*argument) + 320*(r(i,0)*r(i,0)) + 320*(r(i,1)*r(i,1)) + 320*(r(i,2)*r(i,2)) + 4*((r(i,0)*r(i,0)) + (r(i,1)*r(i,1)) + (r(i,2)*r(i,2)))*((alpha*alpha)*rSingleParticle - 20*alpha*argument + 80))*exp(-1.0L/4.0L*alpha*argument)/(argument*argument*argument);
	}

}
