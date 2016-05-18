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
		return 42.1351093655819*exp(-2782.160055*rSingleParticle) + 40.7512433799665*exp(-506.773927*rSingleParticle) + 12.7003825911003*exp(-137.1528019*rSingleParticle);
	case 1:
		return -4.26198581434887*exp(-233.9514118*rSingleParticle) + 5.7007503368583*exp(-54.36527681*rSingleParticle) + 4.30244132054598*exp(-17.68127533*rSingleParticle);
	case 2:
		return -1.61833582599735*exp(-21.45684671*rSingleParticle) + 0.634438826182152*exp(-6.545022156*rSingleParticle) + 1.3087529173072*exp(-2.525273021*rSingleParticle);
	case 3:
		return -0.311678884979035*exp(-1.590049336*rSingleParticle) + 0.00936896129368918*exp(-0.5868282053*rSingleParticle) + 0.292784938811297*exp(-0.2591495227*rSingleParticle);
	case 4:
		return 203.347401124456*r(i,0)*exp(-233.9514118*rSingleParticle) + 127.870302655437*r(i,0)*exp(-54.36527681*rSingleParticle) + 20.2568042861766*r(i,0)*exp(-17.68127533*rSingleParticle);
	case 5:
		return 203.347401124456*r(i,1)*exp(-233.9514118*rSingleParticle) + 127.870302655437*r(i,1)*exp(-54.36527681*rSingleParticle) + 20.2568042861766*r(i,1)*exp(-17.68127533*rSingleParticle);
	case 6:
		return 203.347401124456*r(i,2)*exp(-233.9514118*rSingleParticle) + 127.870302655437*r(i,2)*exp(-54.36527681*rSingleParticle) + 20.2568042861766*r(i,2)*exp(-17.68127533*rSingleParticle);
	case 7:
		return 0.32593817299403*r(i,0)*exp(-21.45684671*rSingleParticle) + 8.62147003303248*r(i,0)*exp(-6.545022156*rSingleParticle) + 2.19912632319612*r(i,0)*exp(-2.525273021*rSingleParticle);
	case 8:
		return 0.32593817299403*r(i,1)*exp(-21.45684671*rSingleParticle) + 8.62147003303248*r(i,1)*exp(-6.545022156*rSingleParticle) + 2.19912632319612*r(i,1)*exp(-2.525273021*rSingleParticle);
	case 9:
		return 0.32593817299403*r(i,2)*exp(-21.45684671*rSingleParticle) + 8.62147003303248*r(i,2)*exp(-6.545022156*rSingleParticle) + 2.19912632319612*r(i,2)*exp(-2.525273021*rSingleParticle);
	case 10:
		return -0.309347835456462*r(i,0)*exp(-1.590049336*rSingleParticle) + 0.418419911098353*r(i,0)*exp(-0.5868282053*rSingleParticle) + 0.144929846938806*r(i,0)*exp(-0.2591495227*rSingleParticle);
	case 11:
		return -0.309347835456462*r(i,1)*exp(-1.590049336*rSingleParticle) + 0.418419911098353*r(i,1)*exp(-0.5868282053*rSingleParticle) + 0.144929846938806*r(i,1)*exp(-0.2591495227*rSingleParticle);
	case 12:
		return -0.309347835456462*r(i,2)*exp(-1.590049336*rSingleParticle) + 0.418419911098353*r(i,2)*exp(-0.5868282053*rSingleParticle) + 0.144929846938806*r(i,2)*exp(-0.2591495227*rSingleParticle);
	case 13:
		return (r(i,0)*r(i,0))*(22.3369316010954*exp(9.070295177*rSingleParticle) + 8.34194591071813*exp(23.982119731*rSingleParticle) + 0.688801507519067*exp(28.001868866*rSingleParticle))*exp(-30.527141887*rSingleParticle);
	case 14:
		return (r(i,1)*r(i,1))*(22.3369316010954*exp(9.070295177*rSingleParticle) + 8.34194591071813*exp(23.982119731*rSingleParticle) + 0.688801507519067*exp(28.001868866*rSingleParticle))*exp(-30.527141887*rSingleParticle);
	case 15:
		return (r(i,2)*r(i,2))*(22.3369316010954*exp(9.070295177*rSingleParticle) + 8.34194591071813*exp(23.982119731*rSingleParticle) + 0.688801507519067*exp(28.001868866*rSingleParticle))*exp(-30.527141887*rSingleParticle);
	case 16:
		return r(i,0)*r(i,1)*(67.0107948032863*exp(9.070295177*rSingleParticle) + 25.0258377321544*exp(23.982119731*rSingleParticle) + 2.0664045225572*exp(28.001868866*rSingleParticle))*exp(-30.527141887*rSingleParticle);
	case 17:
		return r(i,0)*r(i,2)*(67.0107948032863*exp(9.070295177*rSingleParticle) + 25.0258377321544*exp(23.982119731*rSingleParticle) + 2.0664045225572*exp(28.001868866*rSingleParticle))*exp(-30.527141887*rSingleParticle);
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
		retvec(0) = -234453.236379957*r(i,0)*exp(-2782.160055*rSingleParticle) - 41303.3352755967*r(i,0)*exp(-506.773927*rSingleParticle) - 3483.78611514278*r(i,0)*exp(-137.1528019*rSingleParticle);
		retvec(1) = -234453.236379957*r(i,1)*exp(-2782.160055*rSingleParticle) - 41303.3352755967*r(i,1)*exp(-506.773927*rSingleParticle) - 3483.78611514278*r(i,1)*exp(-137.1528019*rSingleParticle);
		retvec(2) = -234453.236379957*r(i,2)*exp(-2782.160055*rSingleParticle) - 41303.3352755967*r(i,2)*exp(-506.773927*rSingleParticle) - 3483.78611514278*r(i,2)*exp(-137.1528019*rSingleParticle);
		return retvec;

	case 1:
		retvec(0) = 1994.19519667698*r(i,0)*exp(-233.9514118*rSingleParticle) - 619.845740176004*r(i,0)*exp(-54.36527681*rSingleParticle) - 152.145299159485*r(i,0)*exp(-17.68127533*rSingleParticle);
		retvec(1) = 1994.19519667698*r(i,1)*exp(-233.9514118*rSingleParticle) - 619.845740176004*r(i,1)*exp(-54.36527681*rSingleParticle) - 152.145299159485*r(i,1)*exp(-17.68127533*rSingleParticle);
		retvec(2) = 1994.19519667698*r(i,2)*exp(-233.9514118*rSingleParticle) - 619.845740176004*r(i,2)*exp(-54.36527681*rSingleParticle) - 152.145299159485*r(i,2)*exp(-17.68127533*rSingleParticle);
		return retvec;

	case 2:
		retvec(0) = 69.4487674874526*r(i,0)*exp(-21.45684671*rSingleParticle) - 8.30483234797764*r(i,0)*exp(-6.545022156*rSingleParticle) - 6.60991686646182*r(i,0)*exp(-2.525273021*rSingleParticle);
		retvec(1) = 69.4487674874526*r(i,1)*exp(-21.45684671*rSingleParticle) - 8.30483234797764*r(i,1)*exp(-6.545022156*rSingleParticle) - 6.60991686646182*r(i,1)*exp(-2.525273021*rSingleParticle);
		retvec(2) = 69.4487674874526*r(i,2)*exp(-21.45684671*rSingleParticle) - 8.30483234797764*r(i,2)*exp(-6.545022156*rSingleParticle) - 6.60991686646182*r(i,2)*exp(-2.525273021*rSingleParticle);
		return retvec;

	case 3:
		retvec(0) = 0.991169608212269*r(i,0)*exp(-1.590049336*rSingleParticle) - 0.0109959414830016*r(i,0)*exp(-0.5868282053*rSingleParticle) - 0.151750154293393*r(i,0)*exp(-0.2591495227*rSingleParticle);
		retvec(1) = 0.991169608212269*r(i,1)*exp(-1.590049336*rSingleParticle) - 0.0109959414830016*r(i,1)*exp(-0.5868282053*rSingleParticle) - 0.151750154293393*r(i,1)*exp(-0.2591495227*rSingleParticle);
		retvec(2) = 0.991169608212269*r(i,2)*exp(-1.590049336*rSingleParticle) - 0.0109959414830016*r(i,2)*exp(-0.5868282053*rSingleParticle) - 0.151750154293393*r(i,2)*exp(-0.2591495227*rSingleParticle);
		return retvec;

	case 4:
		retvec(0) = ((-95146.8231578547*(r(i,0)*r(i,0)) + 203.347401124456)*exp(72.04655214*rSingleParticle) + (-13903.4087992826*(r(i,0)*r(i,0)) + 127.870302655437)*exp(251.63268713*rSingleParticle) + (-716.332267779624*(r(i,0)*r(i,0)) + 20.2568042861766)*exp(288.31668861*rSingleParticle))*exp(-305.99796394*rSingleParticle);
		retvec(1) = -r(i,0)*r(i,1)*(95146.8231578547*exp(72.04655214*rSingleParticle) + 13903.4087992826*exp(251.63268713*rSingleParticle) + 716.332267779624*exp(288.31668861*rSingleParticle))*exp(-305.99796394*rSingleParticle);
		retvec(2) = -r(i,0)*r(i,2)*(95146.8231578547*exp(72.04655214*rSingleParticle) + 13903.4087992826*exp(251.63268713*rSingleParticle) + 716.332267779624*exp(288.31668861*rSingleParticle))*exp(-305.99796394*rSingleParticle);
		return retvec;

	case 5:
		retvec(0) = -r(i,0)*r(i,1)*(95146.8231578547*exp(72.04655214*rSingleParticle) + 13903.4087992826*exp(251.63268713*rSingleParticle) + 716.332267779624*exp(288.31668861*rSingleParticle))*exp(-305.99796394*rSingleParticle);
		retvec(1) = ((-95146.8231578547*(r(i,1)*r(i,1)) + 203.347401124456)*exp(72.04655214*rSingleParticle) + (-13903.4087992826*(r(i,1)*r(i,1)) + 127.870302655437)*exp(251.63268713*rSingleParticle) + (-716.332267779624*(r(i,1)*r(i,1)) + 20.2568042861766)*exp(288.31668861*rSingleParticle))*exp(-305.99796394*rSingleParticle);
		retvec(2) = -r(i,1)*r(i,2)*(95146.8231578547*exp(72.04655214*rSingleParticle) + 13903.4087992826*exp(251.63268713*rSingleParticle) + 716.332267779624*exp(288.31668861*rSingleParticle))*exp(-305.99796394*rSingleParticle);
		return retvec;

	case 6:
		retvec(0) = -r(i,0)*r(i,2)*(95146.8231578547*exp(72.04655214*rSingleParticle) + 13903.4087992826*exp(251.63268713*rSingleParticle) + 716.332267779624*exp(288.31668861*rSingleParticle))*exp(-305.99796394*rSingleParticle);
		retvec(1) = -r(i,1)*r(i,2)*(95146.8231578547*exp(72.04655214*rSingleParticle) + 13903.4087992826*exp(251.63268713*rSingleParticle) + 716.332267779624*exp(288.31668861*rSingleParticle))*exp(-305.99796394*rSingleParticle);
		retvec(2) = ((-95146.8231578547*(r(i,2)*r(i,2)) + 203.347401124456)*exp(72.04655214*rSingleParticle) + (-13903.4087992826*(r(i,2)*r(i,2)) + 127.870302655437)*exp(251.63268713*rSingleParticle) + (-716.332267779624*(r(i,2)*r(i,2)) + 20.2568042861766)*exp(288.31668861*rSingleParticle))*exp(-305.99796394*rSingleParticle);
		return retvec;

	case 7:
		retvec(0) = ((-112.855424766975*(r(i,0)*r(i,0)) + 8.62147003303248)*exp(23.982119731*rSingleParticle) + (-13.9872108297407*(r(i,0)*r(i,0)) + 0.32593817299403)*exp(9.070295177*rSingleParticle) + (-11.1067887474762*(r(i,0)*r(i,0)) + 2.19912632319612)*exp(28.001868866*rSingleParticle))*exp(-30.527141887*rSingleParticle);
		retvec(1) = -r(i,0)*r(i,1)*(13.9872108297407*exp(9.070295177*rSingleParticle) + 112.855424766975*exp(23.982119731*rSingleParticle) + 11.1067887474762*exp(28.001868866*rSingleParticle))*exp(-30.527141887*rSingleParticle);
		retvec(2) = -r(i,0)*r(i,2)*(13.9872108297407*exp(9.070295177*rSingleParticle) + 112.855424766975*exp(23.982119731*rSingleParticle) + 11.1067887474762*exp(28.001868866*rSingleParticle))*exp(-30.527141887*rSingleParticle);
		return retvec;

	case 8:
		retvec(0) = -r(i,0)*r(i,1)*(13.9872108297407*exp(9.070295177*rSingleParticle) + 112.855424766975*exp(23.982119731*rSingleParticle) + 11.1067887474762*exp(28.001868866*rSingleParticle))*exp(-30.527141887*rSingleParticle);
		retvec(1) = ((-112.855424766975*(r(i,1)*r(i,1)) + 8.62147003303248)*exp(23.982119731*rSingleParticle) + (-13.9872108297407*(r(i,1)*r(i,1)) + 0.32593817299403)*exp(9.070295177*rSingleParticle) + (-11.1067887474762*(r(i,1)*r(i,1)) + 2.19912632319612)*exp(28.001868866*rSingleParticle))*exp(-30.527141887*rSingleParticle);
		retvec(2) = -r(i,1)*r(i,2)*(13.9872108297407*exp(9.070295177*rSingleParticle) + 112.855424766975*exp(23.982119731*rSingleParticle) + 11.1067887474762*exp(28.001868866*rSingleParticle))*exp(-30.527141887*rSingleParticle);
		return retvec;

	case 9:
		retvec(0) = -r(i,0)*r(i,2)*(13.9872108297407*exp(9.070295177*rSingleParticle) + 112.855424766975*exp(23.982119731*rSingleParticle) + 11.1067887474762*exp(28.001868866*rSingleParticle))*exp(-30.527141887*rSingleParticle);
		retvec(1) = -r(i,1)*r(i,2)*(13.9872108297407*exp(9.070295177*rSingleParticle) + 112.855424766975*exp(23.982119731*rSingleParticle) + 11.1067887474762*exp(28.001868866*rSingleParticle))*exp(-30.527141887*rSingleParticle);
		retvec(2) = ((-112.855424766975*(r(i,2)*r(i,2)) + 8.62147003303248)*exp(23.982119731*rSingleParticle) + (-13.9872108297407*(r(i,2)*r(i,2)) + 0.32593817299403)*exp(9.070295177*rSingleParticle) + (-11.1067887474762*(r(i,2)*r(i,2)) + 2.19912632319612)*exp(28.001868866*rSingleParticle))*exp(-30.527141887*rSingleParticle);
		return retvec;

	case 10:
		retvec(0) = ((-0.491081210983264*(r(i,0)*r(i,0)) + 0.418419911098353)*exp(1.8491988587*rSingleParticle) + (-0.0751170013183511*(r(i,0)*r(i,0)) + 0.144929846938806)*exp(2.1768775413*rSingleParticle) + (0.98375664072117*(r(i,0)*r(i,0)) - 0.309347835456462)*exp(0.845977728*rSingleParticle))*exp(-2.436027064*rSingleParticle);
		retvec(1) = r(i,0)*r(i,1)*(0.98375664072117*exp(0.845977728*rSingleParticle) - 0.491081210983264*exp(1.8491988587*rSingleParticle) - 0.0751170013183511*exp(2.1768775413*rSingleParticle))*exp(-2.436027064*rSingleParticle);
		retvec(2) = r(i,0)*r(i,2)*(0.98375664072117*exp(0.845977728*rSingleParticle) - 0.491081210983264*exp(1.8491988587*rSingleParticle) - 0.0751170013183511*exp(2.1768775413*rSingleParticle))*exp(-2.436027064*rSingleParticle);
		return retvec;

	case 11:
		retvec(0) = r(i,0)*r(i,1)*(0.98375664072117*exp(0.845977728*rSingleParticle) - 0.491081210983264*exp(1.8491988587*rSingleParticle) - 0.0751170013183511*exp(2.1768775413*rSingleParticle))*exp(-2.436027064*rSingleParticle);
		retvec(1) = ((-0.491081210983264*(r(i,1)*r(i,1)) + 0.418419911098353)*exp(1.8491988587*rSingleParticle) + (-0.0751170013183511*(r(i,1)*r(i,1)) + 0.144929846938806)*exp(2.1768775413*rSingleParticle) + (0.98375664072117*(r(i,1)*r(i,1)) - 0.309347835456462)*exp(0.845977728*rSingleParticle))*exp(-2.436027064*rSingleParticle);
		retvec(2) = r(i,1)*r(i,2)*(0.98375664072117*exp(0.845977728*rSingleParticle) - 0.491081210983264*exp(1.8491988587*rSingleParticle) - 0.0751170013183511*exp(2.1768775413*rSingleParticle))*exp(-2.436027064*rSingleParticle);
		return retvec;

	case 12:
		retvec(0) = r(i,0)*r(i,2)*(0.98375664072117*exp(0.845977728*rSingleParticle) - 0.491081210983264*exp(1.8491988587*rSingleParticle) - 0.0751170013183511*exp(2.1768775413*rSingleParticle))*exp(-2.436027064*rSingleParticle);
		retvec(1) = r(i,1)*r(i,2)*(0.98375664072117*exp(0.845977728*rSingleParticle) - 0.491081210983264*exp(1.8491988587*rSingleParticle) - 0.0751170013183511*exp(2.1768775413*rSingleParticle))*exp(-2.436027064*rSingleParticle);
		retvec(2) = ((-0.491081210983264*(r(i,2)*r(i,2)) + 0.418419911098353)*exp(1.8491988587*rSingleParticle) + (-0.0751170013183511*(r(i,2)*r(i,2)) + 0.144929846938806)*exp(2.1768775413*rSingleParticle) + (0.98375664072117*(r(i,2)*r(i,2)) - 0.309347835456462)*exp(0.845977728*rSingleParticle))*exp(-2.436027064*rSingleParticle);
		return retvec;

	case 13:
		retvec(0) = r(i,0)*((-958.560234672919*(r(i,0)*r(i,0)) + 44.6738632021909)*exp(9.070295177*rSingleParticle) + (-109.196441619607*(r(i,0)*r(i,0)) + 16.6838918214363)*exp(23.982119731*rSingleParticle) + (-3.47882372752406*(r(i,0)*r(i,0)) + 1.37760301503813)*exp(28.001868866*rSingleParticle))*exp(-30.527141887*rSingleParticle);
		retvec(1) = -(r(i,0)*r(i,0))*r(i,1)*(958.560234672919*exp(9.070295177*rSingleParticle) + 109.196441619607*exp(23.982119731*rSingleParticle) + 3.47882372752406*exp(28.001868866*rSingleParticle))*exp(-30.527141887*rSingleParticle);
		retvec(2) = -(r(i,0)*r(i,0))*r(i,2)*(958.560234672919*exp(9.070295177*rSingleParticle) + 109.196441619607*exp(23.982119731*rSingleParticle) + 3.47882372752406*exp(28.001868866*rSingleParticle))*exp(-30.527141887*rSingleParticle);
		return retvec;

	case 14:
		retvec(0) = -r(i,0)*(r(i,1)*r(i,1))*(958.560234672919*exp(9.070295177*rSingleParticle) + 109.196441619607*exp(23.982119731*rSingleParticle) + 3.47882372752406*exp(28.001868866*rSingleParticle))*exp(-30.527141887*rSingleParticle);
		retvec(1) = r(i,1)*((-958.560234672919*(r(i,1)*r(i,1)) + 44.6738632021909)*exp(9.070295177*rSingleParticle) + (-109.196441619607*(r(i,1)*r(i,1)) + 16.6838918214363)*exp(23.982119731*rSingleParticle) + (-3.47882372752406*(r(i,1)*r(i,1)) + 1.37760301503813)*exp(28.001868866*rSingleParticle))*exp(-30.527141887*rSingleParticle);
		retvec(2) = -(r(i,1)*r(i,1))*r(i,2)*(958.560234672919*exp(9.070295177*rSingleParticle) + 109.196441619607*exp(23.982119731*rSingleParticle) + 3.47882372752406*exp(28.001868866*rSingleParticle))*exp(-30.527141887*rSingleParticle);
		return retvec;

	case 15:
		retvec(0) = -r(i,0)*(r(i,2)*r(i,2))*(958.560234672919*exp(9.070295177*rSingleParticle) + 109.196441619607*exp(23.982119731*rSingleParticle) + 3.47882372752406*exp(28.001868866*rSingleParticle))*exp(-30.527141887*rSingleParticle);
		retvec(1) = -r(i,1)*(r(i,2)*r(i,2))*(958.560234672919*exp(9.070295177*rSingleParticle) + 109.196441619607*exp(23.982119731*rSingleParticle) + 3.47882372752406*exp(28.001868866*rSingleParticle))*exp(-30.527141887*rSingleParticle);
		retvec(2) = r(i,2)*((-958.560234672919*(r(i,2)*r(i,2)) + 44.6738632021909)*exp(9.070295177*rSingleParticle) + (-109.196441619607*(r(i,2)*r(i,2)) + 16.6838918214363)*exp(23.982119731*rSingleParticle) + (-3.47882372752406*(r(i,2)*r(i,2)) + 1.37760301503813)*exp(28.001868866*rSingleParticle))*exp(-30.527141887*rSingleParticle);
		return retvec;

	case 16:
		retvec(0) = r(i,1)*((-2875.68070401876*(r(i,0)*r(i,0)) + 67.0107948032863)*exp(9.070295177*rSingleParticle) + (-327.589324858823*(r(i,0)*r(i,0)) + 25.0258377321544)*exp(23.982119731*rSingleParticle) + (-10.4364711825722*(r(i,0)*r(i,0)) + 2.0664045225572)*exp(28.001868866*rSingleParticle))*exp(-30.527141887*rSingleParticle);
		retvec(1) = r(i,0)*((-2875.68070401876*(r(i,1)*r(i,1)) + 67.0107948032863)*exp(9.070295177*rSingleParticle) + (-327.589324858823*(r(i,1)*r(i,1)) + 25.0258377321544)*exp(23.982119731*rSingleParticle) + (-10.4364711825722*(r(i,1)*r(i,1)) + 2.0664045225572)*exp(28.001868866*rSingleParticle))*exp(-30.527141887*rSingleParticle);
		retvec(2) = -r(i,0)*r(i,1)*r(i,2)*(2875.68070401876*exp(9.070295177*rSingleParticle) + 327.589324858823*exp(23.982119731*rSingleParticle) + 10.4364711825722*exp(28.001868866*rSingleParticle))*exp(-30.527141887*rSingleParticle);
		return retvec;

	case 17:
		retvec(0) = r(i,2)*((-2875.68070401876*(r(i,0)*r(i,0)) + 67.0107948032863)*exp(9.070295177*rSingleParticle) + (-327.589324858823*(r(i,0)*r(i,0)) + 25.0258377321544)*exp(23.982119731*rSingleParticle) + (-10.4364711825722*(r(i,0)*r(i,0)) + 2.0664045225572)*exp(28.001868866*rSingleParticle))*exp(-30.527141887*rSingleParticle);
		retvec(1) = -r(i,0)*r(i,1)*r(i,2)*(2875.68070401876*exp(9.070295177*rSingleParticle) + 327.589324858823*exp(23.982119731*rSingleParticle) + 10.4364711825722*exp(28.001868866*rSingleParticle))*exp(-30.527141887*rSingleParticle);
		retvec(2) = r(i,0)*((-2875.68070401876*(r(i,2)*r(i,2)) + 67.0107948032863)*exp(9.070295177*rSingleParticle) + (-327.589324858823*(r(i,2)*r(i,2)) + 25.0258377321544)*exp(23.982119731*rSingleParticle) + (-10.4364711825722*(r(i,2)*r(i,2)) + 2.0664045225572)*exp(28.001868866*rSingleParticle))*exp(-30.527141887*rSingleParticle);
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
		return (955622.053824295*(r(i,0)*r(i,0)) + 955622.053824295*(r(i,1)*r(i,1)) + 955622.053824295*(r(i,2)*r(i,2)) - 10451.3583454283)*exp(-137.1528019*rSingleParticle) + (41862906.8316236*(r(i,0)*r(i,0)) + 41862906.8316236*(r(i,1)*r(i,1)) + 41862906.8316236*(r(i,2)*r(i,2)) - 123910.00582679)*exp(-506.773927*rSingleParticle) + (1304572858.04358*(r(i,0)*r(i,0)) + 1304572858.04358*(r(i,1)*r(i,1)) + 1304572858.04358*(r(i,2)*r(i,2)) - 703359.709139871)*exp(-2782.160055*rSingleParticle);
	case 1:
		return (-933089.563334717*(r(i,0)*r(i,0)) - 933089.563334717*(r(i,1)*r(i,1)) - 933089.563334717*(r(i,2)*r(i,2)) + 5982.58559003094)*exp(-233.9514118*rSingleParticle) + (5380.24584920813*(r(i,0)*r(i,0)) + 5380.24584920813*(r(i,1)*r(i,1)) + 5380.24584920813*(r(i,2)*r(i,2)) - 456.435897478454)*exp(-17.68127533*rSingleParticle) + (67396.1704883356*(r(i,0)*r(i,0)) + 67396.1704883356*(r(i,1)*r(i,1)) + 67396.1704883356*(r(i,2)*r(i,2)) - 1859.53722052801)*exp(-54.36527681*rSingleParticle);
	case 2:
		return (-2980.30311635341*(r(i,0)*r(i,0)) - 2980.30311635341*(r(i,1)*r(i,1)) - 2980.30311635341*(r(i,2)*r(i,2)) + 208.346302462358)*exp(-21.45684671*rSingleParticle) + (33.3836894678578*(r(i,0)*r(i,0)) + 33.3836894678578*(r(i,1)*r(i,1)) + 33.3836894678578*(r(i,2)*r(i,2)) - 19.8297505993855)*exp(-2.525273021*rSingleParticle) + (108.710623438758*(r(i,0)*r(i,0)) + 108.710623438758*(r(i,1)*r(i,1)) + 108.710623438758*(r(i,2)*r(i,2)) - 24.9144970439329)*exp(-6.545022156*rSingleParticle);
	case 3:
		return (-3.1520171548026*(r(i,0)*r(i,0)) - 3.1520171548026*(r(i,1)*r(i,1)) - 3.1520171548026*(r(i,2)*r(i,2)) + 2.97350882463681)*exp(-1.590049336*rSingleParticle) + (0.0129054572121073*(r(i,0)*r(i,0)) + 0.0129054572121073*(r(i,1)*r(i,1)) + 0.0129054572121073*(r(i,2)*r(i,2)) - 0.0329878244490047)*exp(-0.5868282053*rSingleParticle) + (0.0786519601095681*(r(i,0)*r(i,0)) + 0.0786519601095681*(r(i,1)*r(i,1)) + 0.0786519601095681*(r(i,2)*r(i,2)) - 0.455250462880178)*exp(-0.2591495227*rSingleParticle);
	case 4:
		return (25331.3361087496*(r(i,0)*r(i,0)*r(i,0)) + 25331.3361087496*r(i,0)*(r(i,1)*r(i,1)) + 25331.3361087496*r(i,0)*(r(i,2)*r(i,2)) - 3581.66133889812*r(i,0))*exp(-17.68127533*rSingleParticle) + (1511725.33595118*(r(i,0)*r(i,0)*r(i,0)) + 1511725.33595118*r(i,0)*(r(i,1)*r(i,1)) + 1511725.33595118*r(i,0)*(r(i,2)*r(i,2)) - 69517.0439964131*r(i,0))*exp(-54.36527681*rSingleParticle) + (44519467.2121301*(r(i,0)*r(i,0)*r(i,0)) + 44519467.2121301*r(i,0)*(r(i,1)*r(i,1)) + 44519467.2121301*r(i,0)*(r(i,2)*r(i,2)) - 475734.115789274*r(i,0))*exp(-233.9514118*rSingleParticle);
	case 5:
		return (25331.3361087496*(r(i,0)*r(i,0))*r(i,1) + 25331.3361087496*(r(i,1)*r(i,1)*r(i,1)) + 25331.3361087496*r(i,1)*(r(i,2)*r(i,2)) - 3581.66133889812*r(i,1))*exp(-17.68127533*rSingleParticle) + (1511725.33595118*(r(i,0)*r(i,0))*r(i,1) + 1511725.33595118*(r(i,1)*r(i,1)*r(i,1)) + 1511725.33595118*r(i,1)*(r(i,2)*r(i,2)) - 69517.0439964131*r(i,1))*exp(-54.36527681*rSingleParticle) + (44519467.2121301*(r(i,0)*r(i,0))*r(i,1) + 44519467.2121301*(r(i,1)*r(i,1)*r(i,1)) + 44519467.2121301*r(i,1)*(r(i,2)*r(i,2)) - 475734.115789274*r(i,1))*exp(-233.9514118*rSingleParticle);
	case 6:
		return (25331.3361087496*(r(i,0)*r(i,0))*r(i,2) + 25331.3361087496*(r(i,1)*r(i,1))*r(i,2) + 25331.3361087496*(r(i,2)*r(i,2)*r(i,2)) - 3581.66133889812*r(i,2))*exp(-17.68127533*rSingleParticle) + (1511725.33595118*(r(i,0)*r(i,0))*r(i,2) + 1511725.33595118*(r(i,1)*r(i,1))*r(i,2) + 1511725.33595118*(r(i,2)*r(i,2)*r(i,2)) - 69517.0439964131*r(i,2))*exp(-54.36527681*rSingleParticle) + (44519467.2121301*(r(i,0)*r(i,0))*r(i,2) + 44519467.2121301*(r(i,1)*r(i,1))*r(i,2) + 44519467.2121301*(r(i,2)*r(i,2)*r(i,2)) - 475734.115789274*r(i,2))*exp(-233.9514118*rSingleParticle);
	case 7:
		return (56.0953479478959*(r(i,0)*r(i,0)*r(i,0)) + 56.0953479478959*r(i,0)*(r(i,1)*r(i,1)) + 56.0953479478959*r(i,0)*(r(i,2)*r(i,2)) - 55.5339437373809*r(i,0))*exp(-2.525273021*rSingleParticle) + (600.242877348397*(r(i,0)*r(i,0)*r(i,0)) + 600.242877348397*r(i,0)*(r(i,1)*r(i,1)) + 600.242877348397*r(i,0)*(r(i,2)*r(i,2)) - 69.9360541487036*r(i,0))*exp(-21.45684671*rSingleParticle) + (1477.28251104929*(r(i,0)*r(i,0)*r(i,0)) + 1477.28251104929*r(i,0)*(r(i,1)*r(i,1)) + 1477.28251104929*r(i,0)*(r(i,2)*r(i,2)) - 564.277123834876*r(i,0))*exp(-6.545022156*rSingleParticle);
	case 8:
		return (56.0953479478959*(r(i,0)*r(i,0))*r(i,1) + 56.0953479478959*(r(i,1)*r(i,1)*r(i,1)) + 56.0953479478959*r(i,1)*(r(i,2)*r(i,2)) - 55.5339437373809*r(i,1))*exp(-2.525273021*rSingleParticle) + (600.242877348397*(r(i,0)*r(i,0))*r(i,1) + 600.242877348397*(r(i,1)*r(i,1)*r(i,1)) + 600.242877348397*r(i,1)*(r(i,2)*r(i,2)) - 69.9360541487036*r(i,1))*exp(-21.45684671*rSingleParticle) + (1477.28251104929*(r(i,0)*r(i,0))*r(i,1) + 1477.28251104929*(r(i,1)*r(i,1)*r(i,1)) + 1477.28251104929*r(i,1)*(r(i,2)*r(i,2)) - 564.277123834876*r(i,1))*exp(-6.545022156*rSingleParticle);
	case 9:
		return (56.0953479478959*(r(i,0)*r(i,0))*r(i,2) + 56.0953479478959*(r(i,1)*r(i,1))*r(i,2) + 56.0953479478959*(r(i,2)*r(i,2)*r(i,2)) - 55.5339437373808*r(i,2))*exp(-2.525273021*rSingleParticle) + (600.242877348397*(r(i,0)*r(i,0))*r(i,2) + 600.242877348397*(r(i,1)*r(i,1))*r(i,2) + 600.242877348397*(r(i,2)*r(i,2)*r(i,2)) - 69.9360541487036*r(i,2))*exp(-21.45684671*rSingleParticle) + (1477.28251104929*(r(i,0)*r(i,0))*r(i,2) + 1477.28251104929*(r(i,1)*r(i,1))*r(i,2) + 1477.28251104929*(r(i,2)*r(i,2)*r(i,2)) - 564.277123834876*r(i,2))*exp(-6.545022156*rSingleParticle);
	case 10:
		return (-3.12844318672857*(r(i,0)*r(i,0)*r(i,0)) - 3.12844318672857*r(i,0)*(r(i,1)*r(i,1)) - 3.12844318672857*r(i,0)*(r(i,2)*r(i,2)) + 4.91878320360585*r(i,0))*exp(-1.590049336*rSingleParticle) + (0.0389330700766119*(r(i,0)*r(i,0)*r(i,0)) + 0.0389330700766119*r(i,0)*(r(i,1)*r(i,1)) + 0.0389330700766119*r(i,0)*(r(i,2)*r(i,2)) - 0.375585006591756*r(i,0))*exp(-0.2591495227*rSingleParticle) + (0.576360611395719*(r(i,0)*r(i,0)*r(i,0)) + 0.576360611395719*r(i,0)*(r(i,1)*r(i,1)) + 0.576360611395719*r(i,0)*(r(i,2)*r(i,2)) - 2.45540605491632*r(i,0))*exp(-0.5868282053*rSingleParticle);
	case 11:
		return (-3.12844318672857*(r(i,0)*r(i,0))*r(i,1) - 3.12844318672857*(r(i,1)*r(i,1)*r(i,1)) - 3.12844318672857*r(i,1)*(r(i,2)*r(i,2)) + 4.91878320360585*r(i,1))*exp(-1.590049336*rSingleParticle) + (0.0389330700766119*(r(i,0)*r(i,0))*r(i,1) + 0.0389330700766119*(r(i,1)*r(i,1)*r(i,1)) + 0.0389330700766119*r(i,1)*(r(i,2)*r(i,2)) - 0.375585006591756*r(i,1))*exp(-0.2591495227*rSingleParticle) + (0.576360611395719*(r(i,0)*r(i,0))*r(i,1) + 0.576360611395719*(r(i,1)*r(i,1)*r(i,1)) + 0.576360611395719*r(i,1)*(r(i,2)*r(i,2)) - 2.45540605491632*r(i,1))*exp(-0.5868282053*rSingleParticle);
	case 12:
		return (-3.12844318672857*(r(i,0)*r(i,0))*r(i,2) - 3.12844318672857*(r(i,1)*r(i,1))*r(i,2) - 3.12844318672857*(r(i,2)*r(i,2)*r(i,2)) + 4.91878320360585*r(i,2))*exp(-1.590049336*rSingleParticle) + (0.0389330700766119*(r(i,0)*r(i,0))*r(i,2) + 0.0389330700766119*(r(i,1)*r(i,1))*r(i,2) + 0.0389330700766119*(r(i,2)*r(i,2)*r(i,2)) - 0.375585006591756*r(i,2))*exp(-0.2591495227*rSingleParticle) + (0.576360611395719*(r(i,0)*r(i,0))*r(i,2) + 0.576360611395719*(r(i,1)*r(i,1))*r(i,2) + 0.576360611395719*(r(i,2)*r(i,2)*r(i,2)) - 2.45540605491632*r(i,2))*exp(-0.5868282053*rSingleParticle);
	case 13:
		return (17.5699594078623*(r(i,0)*r(i,0)*r(i,0)*r(i,0)) + 17.5699594078623*(r(i,0)*r(i,0))*(r(i,1)*r(i,1)) + 17.5699594078623*(r(i,0)*r(i,0))*(r(i,2)*r(i,2)) - 24.3517660926684*(r(i,0)*r(i,0)) + 1.37760301503813)*exp(-2.525273021*rSingleParticle) + (1429.38625951338*(r(i,0)*r(i,0)*r(i,0)*r(i,0)) + 1429.38625951338*(r(i,0)*r(i,0))*(r(i,1)*r(i,1)) + 1429.38625951338*(r(i,0)*r(i,0))*(r(i,2)*r(i,2)) - 764.375091337252*(r(i,0)*r(i,0)) + 16.6838918214363)*exp(-6.545022156*rSingleParticle) + (41135.3600353569*(r(i,0)*r(i,0)*r(i,0)*r(i,0)) + 41135.3600353569*(r(i,0)*r(i,0))*(r(i,1)*r(i,1)) + 41135.3600353569*(r(i,0)*r(i,0))*(r(i,2)*r(i,2)) - 6709.92164271043*(r(i,0)*r(i,0)) + 44.6738632021909)*exp(-21.45684671*rSingleParticle);
	case 14:
		return (17.5699594078623*(r(i,0)*r(i,0))*(r(i,1)*r(i,1)) + 17.5699594078623*(r(i,1)*r(i,1)*r(i,1)*r(i,1)) + 17.5699594078623*(r(i,1)*r(i,1))*(r(i,2)*r(i,2)) - 24.3517660926684*(r(i,1)*r(i,1)) + 1.37760301503813)*exp(-2.525273021*rSingleParticle) + (1429.38625951338*(r(i,0)*r(i,0))*(r(i,1)*r(i,1)) + 1429.38625951338*(r(i,1)*r(i,1)*r(i,1)*r(i,1)) + 1429.38625951338*(r(i,1)*r(i,1))*(r(i,2)*r(i,2)) - 764.375091337252*(r(i,1)*r(i,1)) + 16.6838918214363)*exp(-6.545022156*rSingleParticle) + (41135.3600353569*(r(i,0)*r(i,0))*(r(i,1)*r(i,1)) + 41135.3600353569*(r(i,1)*r(i,1)*r(i,1)*r(i,1)) + 41135.3600353569*(r(i,1)*r(i,1))*(r(i,2)*r(i,2)) - 6709.92164271043*(r(i,1)*r(i,1)) + 44.6738632021909)*exp(-21.45684671*rSingleParticle);
	case 15:
		return (17.5699594078623*(r(i,0)*r(i,0))*(r(i,2)*r(i,2)) + 17.5699594078623*(r(i,1)*r(i,1))*(r(i,2)*r(i,2)) + 17.5699594078623*(r(i,2)*r(i,2)*r(i,2)*r(i,2)) - 24.3517660926684*(r(i,2)*r(i,2)) + 1.37760301503813)*exp(-2.525273021*rSingleParticle) + (1429.38625951338*(r(i,0)*r(i,0))*(r(i,2)*r(i,2)) + 1429.38625951338*(r(i,1)*r(i,1))*(r(i,2)*r(i,2)) + 1429.38625951338*(r(i,2)*r(i,2)*r(i,2)*r(i,2)) - 764.375091337252*(r(i,2)*r(i,2)) + 16.6838918214363)*exp(-6.545022156*rSingleParticle) + (41135.3600353569*(r(i,0)*r(i,0))*(r(i,2)*r(i,2)) + 41135.3600353569*(r(i,1)*r(i,1))*(r(i,2)*r(i,2)) + 41135.3600353569*(r(i,2)*r(i,2)*r(i,2)*r(i,2)) - 6709.92164271043*(r(i,2)*r(i,2)) + 44.6738632021909)*exp(-21.45684671*rSingleParticle);
	case 16:
		return (52.709878223587*(r(i,0)*r(i,0)*r(i,0))*r(i,1) + 52.709878223587*r(i,0)*(r(i,1)*r(i,1)*r(i,1)) + 52.709878223587*r(i,0)*r(i,1)*(r(i,2)*r(i,2)) - 73.0552982780052*r(i,0)*r(i,1))*exp(-2.525273021*rSingleParticle) + (4288.15877854015*(r(i,0)*r(i,0)*r(i,0))*r(i,1) + 4288.15877854015*r(i,0)*(r(i,1)*r(i,1)*r(i,1)) + 4288.15877854015*r(i,0)*r(i,1)*(r(i,2)*r(i,2)) - 2293.12527401176*r(i,0)*r(i,1))*exp(-6.545022156*rSingleParticle) + (123406.080106071*(r(i,0)*r(i,0)*r(i,0))*r(i,1) + 123406.080106071*r(i,0)*(r(i,1)*r(i,1)*r(i,1)) + 123406.080106071*r(i,0)*r(i,1)*(r(i,2)*r(i,2)) - 20129.7649281313*r(i,0)*r(i,1))*exp(-21.45684671*rSingleParticle);
	case 17:
		return (52.709878223587*(r(i,0)*r(i,0)*r(i,0))*r(i,2) + 52.709878223587*r(i,0)*(r(i,1)*r(i,1))*r(i,2) + 52.709878223587*r(i,0)*(r(i,2)*r(i,2)*r(i,2)) - 73.0552982780052*r(i,0)*r(i,2))*exp(-2.525273021*rSingleParticle) + (4288.15877854015*(r(i,0)*r(i,0)*r(i,0))*r(i,2) + 4288.15877854015*r(i,0)*(r(i,1)*r(i,1))*r(i,2) + 4288.15877854015*r(i,0)*(r(i,2)*r(i,2)*r(i,2)) - 2293.12527401176*r(i,0)*r(i,2))*exp(-6.545022156*rSingleParticle) + (123406.080106071*(r(i,0)*r(i,0)*r(i,0))*r(i,2) + 123406.080106071*r(i,0)*(r(i,1)*r(i,1))*r(i,2) + 123406.080106071*r(i,0)*(r(i,2)*r(i,2)*r(i,2)) - 20129.7649281313*r(i,0)*r(i,2))*exp(-21.45684671*rSingleParticle);
	}

}
