#ifndef MCINTEGRATOR_H
#define MCINTEGRATOR_H

#include <armadillo>

#include "WaveFunction.h"

using namespace arma;

class MCIntegrator
{
public:
    MCIntegrator();

    double runMCIntegration(double a,double b);

private:

    WaveFunction* wf;

    double getGreensFunctionRatio(const mat &y, const mat &x, const mat &quantumForceNew, const mat &quantumForceOld);

    int nDimensions;
    int charge;
    double stepLength;
    int nParticles;


    long idum;

    double alpha;
    double beta;

    int nCycles;
    
    double timestep;
    double D;

    double R;

    mat rOld;
    mat rNew;
};

#endif // MCINTEGRATOR_H
