#include "mcintegrator.h"
#include "LocalEnergy.h"
#include "lib.h"

#include <armadillo>
#include <iostream>

using namespace arma;
using namespace std;


MCIntegrator:: MCIntegrator() :
    nDimensions(3),
    charge(4),
    stepLength(0.1),
    nParticles(4),
    h(0.001),
    h2(1000000),
    idum(-1),
    alpha(3.3),
    beta(0.3),
    nCycles(1000000)
{
}

void MCIntegrator::runMCIntegration() {

  LocalEnergy* le = new LocalEnergy();

  int accetpted = 0;

  double waveFunctionOld = 0;
  double waveFunctionNew = 0;

  double energySum = 0;
  double energySquaredSum = 0;

  double deltaE;

  rOld = zeros<mat>( nParticles , nDimensions );
  rNew = zeros<mat>( nParticles , nDimensions );
 
  for ( int i = 0; i < nParticles ; i ++) {
    for ( int j = 0; j < nDimensions ; j ++) {
      rOld (i , j) = stepLength * ( ran2 (&idum) - 0.5) ;
    }
  } 

 
  waveFunctionOld = slater(rOld);

  for(int cycle = 0; cycle < nCycles; cycle++) {

    // New position to test
    for(int i = 0; i < nParticles; i++) {
      for(int j = 0; j < nDimensions; j++) {
	rNew(i,j) = rOld(i,j) + stepLength*(ran2(&idum) - 0.5);
      }
    }

    // Recalculate the value of the wave function
    waveFunctionNew = slater(rNew);
    
    // Metropolis test
    if (ran2(&idum) <= (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld)) {

      rOld = rNew;

      waveFunctionOld = waveFunctionNew;
	
      accetpted++;
	
    }

    // update energies
    deltaE = localEnergy(rOld);
    
    //deltaE = le->get_local_energy(rNew,nDimensions,nParticles,charge,alpha,beta);
    
    
    energySum += deltaE;
    energySquaredSum += deltaE*deltaE;
  }

  double energy = energySum/(nCycles);// * nParticles);
  double energySquared = energySquaredSum/(nCycles);// * nParticles);
  cout << "Energy: " << energy << " Energy (squared sum): " << energySquared << endl;
  cout << "Acceptred steps: " << (double)accetpted/(double)(nCycles * nParticles) << endl;        

}

double MCIntegrator::slater(const mat &r) {

  mat detPlus  = zeros<mat>( 2 , 2 );
  mat detMinus = zeros<mat>( 2 , 2 );

  detPlus(0,0) = waveFunction1s(r,0);
  detPlus(1,0) = waveFunction1s(r,1);
  detPlus(0,1) = waveFunction2s(r,0);
  detPlus(1,1) = waveFunction2s(r,1);

  detMinus(0,0) = waveFunction1s(r,2);
  detMinus(1,0) = waveFunction1s(r,3);
  detMinus(0,1) = waveFunction2s(r,2);
  detMinus(1,1) = waveFunction2s(r,3);

  
  double jast = 1.0;
  
  for(int i = 0; i < nParticles; i++) {
    for(int j = 0; j < i; j++) {
      jast *= jastrowFactor(r,i,j);
    }
  }
 
  return det(detPlus)*det(detMinus) * jast;

}

double MCIntegrator::waveFunction1s(const mat &r,int nParticle) {
  
  int i = nParticle;
  double rSingleParticle = 0;
  double argument = 0.0;
  
  for(int j = 0; j < nDimensions; j++) {
    rSingleParticle += r(i,j) * r(i,j);
  }
  argument += sqrt(rSingleParticle);
  return exp(-argument * alpha);// * jastrowFactor(r);
  
}

double MCIntegrator::waveFunction2s(const mat &r,int nParticle) {
  
  int i = nParticle;
  double rSingleParticle = 0;
  double argument = 0.0;
  
  for(int j = 0; j < nDimensions; j++) {
    rSingleParticle += r(i,j) * r(i,j);
  }
  argument += sqrt(rSingleParticle);
  return (1.0 - alpha*argument/2.0)*exp(-argument * alpha/2.0);// * jastrowFactor(r);
  
}

double MCIntegrator::jastrowFactor(const mat &r,int i,int j) {
  
  rowvec r12 = r.row(i) - r.row(j);
  double r12norm = norm(r12, 2);
  
  return exp(alpha*r12norm / (2 * (1 + beta * r12norm)));
}

/*
double MCIntegrator::localEnergyAnalytical(const mat &r) {

  double r1 = norm(r.row(0),2);//sqrt(lngd1);
  double r2 = norm(r.row(1),2);//sqrt(lngd2);

  double r12 = norm(r.row(1) - r.row(0),2);

  double dott = dot(r.row(0), r.row(1));  

  double div = 1 + beta*r12;
  double div2 = div*div;

  double EL1 = (alpha - charge)*(1/r1 + 1/r2) + 1/r12 - alpha*alpha;

  double EL = (alpha*(r1 + r2)/r12*(1-dott/(r1*r2)) - 1/(2*div2) - 2/r12 + 2*beta/div)/(2*div2);

  return EL1 + EL;

}
*/
double MCIntegrator::localEnergy(const mat &r) {

    mat rPlus = zeros<mat>(nParticles, nDimensions);
    mat rMinus = zeros<mat>(nParticles, nDimensions);

    rPlus = rMinus = r;

    double waveFunctionMinus = 0;
    double waveFunctionPlus = 0;

    double waveFunctionCurrent = slater(r);

    // Kinetic energy

    double kineticEnergy = 0;
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rPlus(i,j) += h;
            rMinus(i,j) -= h;
            waveFunctionMinus = slater(rMinus);
            waveFunctionPlus  = slater(rPlus);
            kineticEnergy -= (waveFunctionMinus + waveFunctionPlus - 2 * waveFunctionCurrent);
            rPlus(i,j) = r(i,j);
            rMinus(i,j) = r(i,j);
        }
    }
    kineticEnergy = 0.5 * h2 * kineticEnergy / waveFunctionCurrent;

    // Potential energy
    double potentialEnergy = 0;
    double rSingleParticle = 0;
    for(int i = 0; i < nParticles; i++) {
        rSingleParticle = 0;
        for(int j = 0; j < nDimensions; j++) {
            rSingleParticle += r(i,j)*r(i,j);
        }
        potentialEnergy -= charge / sqrt(rSingleParticle);
    }
    // Contribution from electron-electron potential
    double r12 = 0;
    for(int i = 0; i < nParticles; i++) {
        for(int j = i + 1; j < nParticles; j++) {
            r12 = 0;
            for(int k = 0; k < nDimensions; k++) {
                r12 += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
            }
	    potentialEnergy += 1 / sqrt(r12);
        }
    }

    return kineticEnergy + potentialEnergy;
}
