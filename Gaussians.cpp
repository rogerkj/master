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
		return 4.25194331812225*exp(-130.70932*rSingleParticle) + 4.11229375245432*exp(-23.808861*rSingleParticle) + 1.28162254318487*exp(-6.44360829999999*rSingleParticle);
	case 1:
		return -0.23941300621163*exp(-5.0331513*rSingleParticle) + 0.320234233437093*exp(-1.1695961*rSingleParticle) + 0.241685574000907*exp(-0.380389*rSingleParticle);
	case 2:
		return 1.67545013446767*r(i,0)*exp(-5.0331513*rSingleParticle) + 1.05356801827835*r(i,0)*exp(-1.1695961*rSingleParticle) + 0.166902899704829*r(i,0)*exp(-0.380389*rSingleParticle);
	case 3:
		return 1.67545013446767*r(i,1)*exp(-5.0331513*rSingleParticle) + 1.05356801827835*r(i,1)*exp(-1.1695961*rSingleParticle) + 0.166902899704829*r(i,1)*exp(-0.380389*rSingleParticle);
	case 4:
		return 1.67545013446767*r(i,2)*exp(-5.0331513*rSingleParticle) + 1.05356801827835*r(i,2)*exp(-1.1695961*rSingleParticle) + 0.166902899704829*r(i,2)*exp(-0.380389*rSingleParticle);
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
		retvec(0) = -1111.53723958061*r(i,0)*exp(-130.70932*rSingleParticle) - 195.818060686707*r(i,0)*exp(-23.808861*rSingleParticle) - 16.5165473134663*r(i,0)*exp(-6.44360829999999*rSingleParticle);
		retvec(1) = -1111.53723958061*r(i,1)*exp(-130.70932*rSingleParticle) - 195.818060686707*r(i,1)*exp(-23.808861*rSingleParticle) - 16.5165473134663*r(i,1)*exp(-6.44360829999999*rSingleParticle);
		retvec(2) = -1111.53723958061*r(i,2)*exp(-130.70932*rSingleParticle) - 195.818060686707*r(i,2)*exp(-23.808861*rSingleParticle) - 16.5165473134663*r(i,2)*exp(-6.44360829999999*rSingleParticle);
		return retvec;

	case 1:
		retvec(0) = 2.41000376690194*r(i,0)*exp(-5.0331513*rSingleParticle) - 0.749089421029026*r(i,0)*exp(-1.1695961*rSingleParticle) - 0.183869067617262*r(i,0)*exp(-0.380389*rSingleParticle);
		retvec(1) = 2.41000376690194*r(i,1)*exp(-5.0331513*rSingleParticle) - 0.749089421029026*r(i,1)*exp(-1.1695961*rSingleParticle) - 0.183869067617262*r(i,1)*exp(-0.380389*rSingleParticle);
		retvec(2) = 2.41000376690194*r(i,2)*exp(-5.0331513*rSingleParticle) - 0.749089421029026*r(i,2)*exp(-1.1695961*rSingleParticle) - 0.183869067617262*r(i,2)*exp(-0.380389*rSingleParticle);
		return retvec;

	case 2:
		retvec(0) = ((-16.8655880447622*(r(i,0)*r(i,0)) + 1.67545013446767)*exp(1.5499851*rSingleParticle) + (-2.46449809052618*(r(i,0)*r(i,0)) + 1.05356801827835)*exp(5.4135403*rSingleParticle) + (-0.126976054231641*(r(i,0)*r(i,0)) + 0.166902899704829)*exp(6.2027474*rSingleParticle))*exp(-6.5831364*rSingleParticle);
		retvec(1) = -r(i,0)*r(i,1)*(16.8655880447622*exp(1.5499851*rSingleParticle) + 2.46449809052618*exp(5.4135403*rSingleParticle) + 0.126976054231641*exp(6.2027474*rSingleParticle))*exp(-6.5831364*rSingleParticle);
		retvec(2) = -r(i,0)*r(i,2)*(16.8655880447622*exp(1.5499851*rSingleParticle) + 2.46449809052618*exp(5.4135403*rSingleParticle) + 0.126976054231641*exp(6.2027474*rSingleParticle))*exp(-6.5831364*rSingleParticle);
		return retvec;

	case 3:
		retvec(0) = -r(i,0)*r(i,1)*(16.8655880447622*exp(1.5499851*rSingleParticle) + 2.46449809052618*exp(5.4135403*rSingleParticle) + 0.126976054231641*exp(6.2027474*rSingleParticle))*exp(-6.5831364*rSingleParticle);
		retvec(1) = ((-16.8655880447622*(r(i,1)*r(i,1)) + 1.67545013446767)*exp(1.5499851*rSingleParticle) + (-2.46449809052618*(r(i,1)*r(i,1)) + 1.05356801827835)*exp(5.4135403*rSingleParticle) + (-0.126976054231641*(r(i,1)*r(i,1)) + 0.166902899704829)*exp(6.2027474*rSingleParticle))*exp(-6.5831364*rSingleParticle);
		retvec(2) = -r(i,1)*r(i,2)*(16.8655880447622*exp(1.5499851*rSingleParticle) + 2.46449809052618*exp(5.4135403*rSingleParticle) + 0.126976054231641*exp(6.2027474*rSingleParticle))*exp(-6.5831364*rSingleParticle);
		return retvec;

	case 4:
		retvec(0) = -r(i,0)*r(i,2)*(16.8655880447622*exp(1.5499851*rSingleParticle) + 2.46449809052618*exp(5.4135403*rSingleParticle) + 0.126976054231641*exp(6.2027474*rSingleParticle))*exp(-6.5831364*rSingleParticle);
		retvec(1) = -r(i,1)*r(i,2)*(16.8655880447622*exp(1.5499851*rSingleParticle) + 2.46449809052618*exp(5.4135403*rSingleParticle) + 0.126976054231641*exp(6.2027474*rSingleParticle))*exp(-6.5831364*rSingleParticle);
		retvec(2) = ((-16.8655880447622*(r(i,2)*r(i,2)) + 1.67545013446767)*exp(1.5499851*rSingleParticle) + (-2.46449809052618*(r(i,2)*r(i,2)) + 1.05356801827835)*exp(5.4135403*rSingleParticle) + (-0.126976054231641*(r(i,2)*r(i,2)) + 0.166902899704829)*exp(6.2027474*rSingleParticle))*exp(-6.5831364*rSingleParticle);
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
		return (212.852322712788*(r(i,0)*r(i,0)) + 212.852322712788*(r(i,1)*r(i,1)) + 212.852322712788*(r(i,2)*r(i,2)) - 49.5496419403988)*exp(-6.4436083*rSingleParticle) + (9324.40997635873*(r(i,0)*r(i,0)) + 9324.40997635873*(r(i,1)*r(i,1)) + 9324.40997635873*(r(i,2)*r(i,2)) - 587.45418206012)*exp(-23.808861*rSingleParticle) + (290576.553480516*(r(i,0)*r(i,0)) + 290576.553480516*(r(i,1)*r(i,1)) + 290576.553480516*(r(i,2)*r(i,2)) - 3334.61171874182)*exp(-130.70932*rSingleParticle);
	case 1:
		return (-24.2598271847748*(r(i,0)*r(i,0)) - 24.2598271847748*(r(i,1)*r(i,1)) - 24.2598271847748*(r(i,2)*r(i,2)) + 7.23001130070583)*exp(-5.0331513*rSingleParticle) + (0.139883541523725*(r(i,0)*r(i,0)) + 0.139883541523725*(r(i,1)*r(i,1)) + 0.139883541523725*(r(i,2)*r(i,2)) - 0.551607202851786)*exp(-0.380389*rSingleParticle) + (1.75226413077361*(r(i,0)*r(i,0)) + 1.75226413077361*(r(i,1)*r(i,1)) + 1.75226413077361*(r(i,2)*r(i,2)) - 2.24726826308708)*exp(-1.1695961*rSingleParticle);
	case 2:
		return (0.096600588586239*(r(i,0)*r(i,0)*r(i,0)) + 0.096600588586239*r(i,0)*(r(i,1)*r(i,1)) + 0.096600588586239*r(i,0)*(r(i,2)*r(i,2)) - 0.634880271158203*r(i,0))*exp(-0.380389*rSingleParticle) + (5.76493471027373*(r(i,0)*r(i,0)*r(i,0)) + 5.76493471027373*r(i,0)*(r(i,1)*r(i,1)) + 5.76493471027373*r(i,0)*(r(i,2)*r(i,2)) - 12.3224904526309*r(i,0))*exp(-1.1695961*rSingleParticle) + (169.774112785519*(r(i,0)*r(i,0)*r(i,0)) + 169.774112785519*r(i,0)*(r(i,1)*r(i,1)) + 169.774112785519*r(i,0)*(r(i,2)*r(i,2)) - 84.3279402238112*r(i,0))*exp(-5.0331513*rSingleParticle);
	case 3:
		return (0.096600588586239*(r(i,0)*r(i,0))*r(i,1) + 0.096600588586239*(r(i,1)*r(i,1)*r(i,1)) + 0.096600588586239*r(i,1)*(r(i,2)*r(i,2)) - 0.634880271158203*r(i,1))*exp(-0.380389*rSingleParticle) + (5.76493471027373*(r(i,0)*r(i,0))*r(i,1) + 5.76493471027373*(r(i,1)*r(i,1)*r(i,1)) + 5.76493471027373*r(i,1)*(r(i,2)*r(i,2)) - 12.3224904526309*r(i,1))*exp(-1.1695961*rSingleParticle) + (169.774112785519*(r(i,0)*r(i,0))*r(i,1) + 169.774112785519*(r(i,1)*r(i,1)*r(i,1)) + 169.774112785519*r(i,1)*(r(i,2)*r(i,2)) - 84.3279402238112*r(i,1))*exp(-5.0331513*rSingleParticle);
	case 4:
		return (0.096600588586239*(r(i,0)*r(i,0))*r(i,2) + 0.096600588586239*(r(i,1)*r(i,1))*r(i,2) + 0.096600588586239*(r(i,2)*r(i,2)*r(i,2)) - 0.634880271158203*r(i,2))*exp(-0.380389*rSingleParticle) + (5.76493471027373*(r(i,0)*r(i,0))*r(i,2) + 5.76493471027373*(r(i,1)*r(i,1))*r(i,2) + 5.76493471027373*(r(i,2)*r(i,2)*r(i,2)) - 12.3224904526309*r(i,2))*exp(-1.1695961*rSingleParticle) + (169.774112785519*(r(i,0)*r(i,0))*r(i,2) + 169.774112785519*(r(i,1)*r(i,1))*r(i,2) + 169.774112785519*(r(i,2)*r(i,2)*r(i,2)) - 84.3279402238112*r(i,2))*exp(-5.0331513*rSingleParticle);
	}

}
