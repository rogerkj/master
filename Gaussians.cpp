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
		return 6.00288289131012*exp(-207.01561*rSingleParticle) + 5.80572677377551*exp(-37.708151*rSingleParticle) + 1.80939141292887*exp(-10.205297*rSingleParticle);
	case 1:
		return -0.346707067600285*exp(-8.2463151*rSingleParticle) + 0.463748697249596*exp(-1.9162662*rSingleParticle) + 0.349998081143196*exp(-0.6232293*rSingleParticle);
	case 2:
		return 3.10567804189332*r(i,0)*exp(-8.2463151*rSingleParticle) + 1.95293365076282*r(i,0)*exp(-1.9162662*rSingleParticle) + 0.309377533950825*r(i,0)*exp(-0.6232293*rSingleParticle);
	case 3:
		return 3.10567804189332*r(i,1)*exp(-8.2463151*rSingleParticle) + 1.95293365076282*r(i,1)*exp(-1.9162662*rSingleParticle) + 0.309377533950825*r(i,1)*exp(-0.6232293*rSingleParticle);
	case 4:
		return 3.10567804189332*r(i,2)*exp(-8.2463151*rSingleParticle) + 1.95293365076282*r(i,2)*exp(-1.9162662*rSingleParticle) + 0.309377533950825*r(i,2)*exp(-0.6232293*rSingleParticle);
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
		retvec(0) = -2485.38092700626*r(i,0)*exp(-207.01561*rSingleParticle) - 437.846443700539*r(i,0)*exp(-37.708151*rSingleParticle) - 36.9307535163774*r(i,0)*exp(-10.205297*rSingleParticle);
		retvec(1) = -2485.38092700626*r(i,1)*exp(-207.01561*rSingleParticle) - 437.846443700539*r(i,1)*exp(-37.708151*rSingleParticle) - 36.9307535163774*r(i,1)*exp(-10.205297*rSingleParticle);
		retvec(2) = -2485.38092700626*r(i,2)*exp(-207.01561*rSingleParticle) - 437.846443700539*r(i,2)*exp(-37.708151*rSingleParticle) - 36.9307535163774*r(i,2)*exp(-10.205297*rSingleParticle);
		return retvec;

	case 1:
		retvec(0) = 5.7181114536579*r(i,0)*exp(-8.2463151*rSingleParticle) - 1.77733190766687*r(i,0)*exp(-1.9162662*rSingleParticle) - 0.436258118224435*r(i,0)*exp(-0.6232293*rSingleParticle);
		retvec(1) = 5.7181114536579*r(i,1)*exp(-8.2463151*rSingleParticle) - 1.77733190766687*r(i,1)*exp(-1.9162662*rSingleParticle) - 0.436258118224435*r(i,1)*exp(-0.6232293*rSingleParticle);
		retvec(2) = 5.7181114536579*r(i,2)*exp(-8.2463151*rSingleParticle) - 1.77733190766687*r(i,2)*exp(-1.9162662*rSingleParticle) - 0.436258118224435*r(i,2)*exp(-0.6232293*rSingleParticle);
		return retvec;

	case 2:
		retvec(0) = ((-51.2207994652066*(r(i,0)*r(i,0)) + 3.10567804189332)*exp(2.5394955*rSingleParticle) + (-7.48468149159878*(r(i,0)*r(i,0)) + 1.95293365076282)*exp(8.8695444*rSingleParticle) + (-0.385626287839798*(r(i,0)*r(i,0)) + 0.309377533950825)*exp(10.1625813*rSingleParticle))*exp(-10.7858106*rSingleParticle);
		retvec(1) = -r(i,0)*r(i,1)*(51.2207994652066*exp(2.5394955*rSingleParticle) + 7.48468149159878*exp(8.8695444*rSingleParticle) + 0.385626287839798*exp(10.1625813*rSingleParticle))*exp(-10.7858106*rSingleParticle);
		retvec(2) = -r(i,0)*r(i,2)*(51.2207994652066*exp(2.5394955*rSingleParticle) + 7.48468149159878*exp(8.8695444*rSingleParticle) + 0.385626287839798*exp(10.1625813*rSingleParticle))*exp(-10.7858106*rSingleParticle);
		return retvec;

	case 3:
		retvec(0) = -r(i,0)*r(i,1)*(51.2207994652066*exp(2.5394955*rSingleParticle) + 7.48468149159878*exp(8.8695444*rSingleParticle) + 0.385626287839798*exp(10.1625813*rSingleParticle))*exp(-10.7858106*rSingleParticle);
		retvec(1) = ((-51.2207994652066*(r(i,1)*r(i,1)) + 3.10567804189332)*exp(2.5394955*rSingleParticle) + (-7.48468149159878*(r(i,1)*r(i,1)) + 1.95293365076282)*exp(8.8695444*rSingleParticle) + (-0.385626287839798*(r(i,1)*r(i,1)) + 0.309377533950825)*exp(10.1625813*rSingleParticle))*exp(-10.7858106*rSingleParticle);
		retvec(2) = -r(i,1)*r(i,2)*(51.2207994652066*exp(2.5394955*rSingleParticle) + 7.48468149159878*exp(8.8695444*rSingleParticle) + 0.385626287839798*exp(10.1625813*rSingleParticle))*exp(-10.7858106*rSingleParticle);
		return retvec;

	case 4:
		retvec(0) = -r(i,0)*r(i,2)*(51.2207994652066*exp(2.5394955*rSingleParticle) + 7.48468149159878*exp(8.8695444*rSingleParticle) + 0.385626287839798*exp(10.1625813*rSingleParticle))*exp(-10.7858106*rSingleParticle);
		retvec(1) = -r(i,1)*r(i,2)*(51.2207994652066*exp(2.5394955*rSingleParticle) + 7.48468149159878*exp(8.8695444*rSingleParticle) + 0.385626287839798*exp(10.1625813*rSingleParticle))*exp(-10.7858106*rSingleParticle);
		retvec(2) = ((-51.2207994652066*(r(i,2)*r(i,2)) + 3.10567804189332)*exp(2.5394955*rSingleParticle) + (-7.48468149159878*(r(i,2)*r(i,2)) + 1.95293365076282)*exp(8.8695444*rSingleParticle) + (-0.385626287839798*(r(i,2)*r(i,2)) + 0.309377533950825)*exp(10.1625813*rSingleParticle))*exp(-10.7858106*rSingleParticle);
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
		return (753.778616136852*(r(i,0)*r(i,0)) + 753.778616136852*(r(i,1)*r(i,1)) + 753.778616136852*(r(i,2)*r(i,2)) - 110.792260549132)*exp(-10.205297*rSingleParticle) + (33020.7596277459*(r(i,0)*r(i,0)) + 33020.7596277459*(r(i,1)*r(i,1)) + 33020.7596277459*(r(i,2)*r(i,2)) - 1313.53933110162)*exp(-37.708151*rSingleParticle) + (1029025.29737313*(r(i,0)*r(i,0)) + 1029025.29737313*(r(i,1)*r(i,1)) + 1029025.29737313*(r(i,2)*r(i,2)) - 7456.14278101877)*exp(-207.01561*rSingleParticle);
	case 1:
		return (-94.3066976475641*(r(i,0)*r(i,0)) - 94.3066976475641*(r(i,1)*r(i,1)) - 94.3066976475641*(r(i,2)*r(i,2)) + 17.1543343609737)*exp(-8.2463151*rSingleParticle) + (0.543777683280663*(r(i,0)*r(i,0)) + 0.543777683280663*(r(i,1)*r(i,1)) + 0.543777683280663*(r(i,2)*r(i,2)) - 1.3087743546733)*exp(-0.6232293*rSingleParticle) + (6.81168212168708*(r(i,0)*r(i,0)) + 6.81168212168708*(r(i,1)*r(i,1)) + 6.81168212168708*(r(i,2)*r(i,2)) - 5.33199572300061)*exp(-1.9162662*rSingleParticle);
	case 2:
		return (0.480667202863991*(r(i,0)*r(i,0)*r(i,0)) + 0.480667202863991*r(i,0)*(r(i,1)*r(i,1)) + 0.480667202863991*r(i,0)*(r(i,2)*r(i,2)) - 1.92813143919899*r(i,0))*exp(-0.6232293*rSingleParticle) + (28.6852843202326*(r(i,0)*r(i,0)*r(i,0)) + 28.6852843202326*r(i,0)*(r(i,1)*r(i,1)) + 28.6852843202326*r(i,0)*(r(i,2)*r(i,2)) - 37.4234074579939*r(i,0))*exp(-1.9162662*rSingleParticle) + (844.76570412801*(r(i,0)*r(i,0)*r(i,0)) + 844.76570412801*r(i,0)*(r(i,1)*r(i,1)) + 844.76570412801*r(i,0)*(r(i,2)*r(i,2)) - 256.103997326033*r(i,0))*exp(-8.2463151*rSingleParticle);
	case 3:
		return (0.480667202863991*(r(i,0)*r(i,0))*r(i,1) + 0.480667202863991*(r(i,1)*r(i,1)*r(i,1)) + 0.480667202863991*r(i,1)*(r(i,2)*r(i,2)) - 1.92813143919899*r(i,1))*exp(-0.6232293*rSingleParticle) + (28.6852843202326*(r(i,0)*r(i,0))*r(i,1) + 28.6852843202326*(r(i,1)*r(i,1)*r(i,1)) + 28.6852843202326*r(i,1)*(r(i,2)*r(i,2)) - 37.4234074579939*r(i,1))*exp(-1.9162662*rSingleParticle) + (844.76570412801*(r(i,0)*r(i,0))*r(i,1) + 844.76570412801*(r(i,1)*r(i,1)*r(i,1)) + 844.76570412801*r(i,1)*(r(i,2)*r(i,2)) - 256.103997326033*r(i,1))*exp(-8.2463151*rSingleParticle);
	case 4:
		return (0.480667202863991*(r(i,0)*r(i,0))*r(i,2) + 0.480667202863991*(r(i,1)*r(i,1))*r(i,2) + 0.480667202863991*(r(i,2)*r(i,2)*r(i,2)) - 1.92813143919899*r(i,2))*exp(-0.6232293*rSingleParticle) + (28.6852843202326*(r(i,0)*r(i,0))*r(i,2) + 28.6852843202326*(r(i,1)*r(i,1))*r(i,2) + 28.6852843202326*(r(i,2)*r(i,2)*r(i,2)) - 37.4234074579939*r(i,2))*exp(-1.9162662*rSingleParticle) + (844.76570412801*(r(i,0)*r(i,0))*r(i,2) + 844.76570412801*(r(i,1)*r(i,1))*r(i,2) + 844.76570412801*(r(i,2)*r(i,2)*r(i,2)) - 256.103997326033*r(i,2))*exp(-8.2463151*rSingleParticle);
	}

}
