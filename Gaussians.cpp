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

	argument += sqrt(rSingleParticle);
	switch (orbital) {
	case 0:
		return 1.4158460928168*exp(-30.167871*rSingleParticle) + 1.36934465925061*exp(-5.4951153*rSingleParticle) + 0.426764928968415*exp(-1.4871927*rSingleParticle);
	case 1:
		return -0.0874824066581176*exp(-1.3148331*rSingleParticle) + 0.117014776148572*exp(-0.3055389*rSingleParticle) + 0.088312774310746*exp(-0.0993706999999999*rSingleParticle);
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
		retvec(0) = -85.4261245679023*r(i,0)*exp(-30.167871*rSingleParticle) - 15.0494135760426*r(i,0)*exp(-5.4951153*rSingleParticle) - 1.26936337395569*r(i,0)*exp(-1.4871927*rSingleParticle);
		retvec(1) = -85.4261245679023*r(i,1)*exp(-30.167871*rSingleParticle) - 15.0494135760426*r(i,1)*exp(-5.4951153*rSingleParticle) - 1.26936337395569*r(i,1)*exp(-1.4871927*rSingleParticle);
		retvec(2) = -85.4261245679023*r(i,2)*exp(-30.167871*rSingleParticle) - 15.0494135760426*r(i,2)*exp(-5.4951153*rSingleParticle) - 1.26936337395569*r(i,2)*exp(-1.4871927*rSingleParticle);
		return retvec;

	case 1:
		retvec(0) = 0.230049527883507*r(i,0)*exp(-1.3148331*rSingleParticle) - 0.0715051319763618*r(i,0)*exp(-0.3055389*rSingleParticle) - 0.0175514044044017*r(i,0)*exp(-0.0993706999999999*rSingleParticle);
		retvec(1) = 0.230049527883507*r(i,1)*exp(-1.3148331*rSingleParticle) - 0.0715051319763618*r(i,1)*exp(-0.3055389*rSingleParticle) - 0.0175514044044017*r(i,1)*exp(-0.0993706999999999*rSingleParticle);
		retvec(2) = 0.230049527883507*r(i,2)*exp(-1.3148331*rSingleParticle) - 0.0715051319763618*r(i,2)*exp(-0.3055389*rSingleParticle) - 0.0175514044044017*r(i,2)*exp(-0.0993706999999999*rSingleParticle);
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
		return (3.77557588678855*(r(i,0)*r(i,0)) + 3.77557588678855*(r(i,1)*r(i,1)) + 3.77557588678855*(r(i,2)*r(i,2)) - 3.80809012186707)*exp(-1.4871927*rSingleParticle) + (165.396525595479*(r(i,0)*r(i,0)) + 165.396525595479*(r(i,1)*r(i,1)) + 165.396525595479*(r(i,2)*r(i,2)) - 45.1482407281277)*exp(-5.4951153*rSingleParticle) + (5154.24861198881*(r(i,0)*r(i,0)) + 5154.24861198881*(r(i,1)*r(i,1)) + 5154.24861198881*(r(i,2)*r(i,2)) - 256.278373703707)*exp(-30.167871*rSingleParticle);
	case 1:
		return (-0.604953467801215*(r(i,0)*r(i,0)) - 0.604953467801215*(r(i,1)*r(i,1)) - 0.604953467801215*(r(i,2)*r(i,2)) + 0.69014858365052)*exp(-1.3148331*rSingleParticle) + (0.00348819068329696*(r(i,0)*r(i,0)) + 0.00348819068329696*(r(i,1)*r(i,1)) + 0.00348819068329696*(r(i,2)*r(i,2)) - 0.0526542132132051)*exp(-0.0993707*rSingleParticle) + (0.0436951987368248*(r(i,0)*r(i,0)) + 0.0436951987368248*(r(i,1)*r(i,1)) + 0.0436951987368248*(r(i,2)*r(i,2)) - 0.214515395929085)*exp(-0.3055389*rSingleParticle);
	}

}
