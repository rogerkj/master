
#include "setup.h"

#include "mcintegrator.h"
#include "LocalEnergyNumeric.h"
#include "LocalEnergyOpt.h"
#include "LocalEnergy.h"
#include "QuantumForce.h"
#include "QuantumForceNumeric.h"
#include "Hydrogenic.h"
#include "Gaussians.h"
#include "Diatomic.h"
#include "lib.h"

#include <armadillo>
#include <iostream>

using namespace arma;
using namespace std;


MCIntegrator:: MCIntegrator() :
    nDimensions(3),
    charge(CHARGE),
    stepLength(0.1),
    nParticles(NUMPART),
    idum(-1),
    nCycles(NUMCYCLES),
    timestep(0.01),
    D(0.5),
    R(4.63)
   
{
}


double MCIntegrator::runMCIntegration(double a,double b) {

  alpha = a;
  beta = b;

#ifndef DIATOMIC

#ifndef GAUSSIAN 

  Orbitals* dr = new Hydrogenic(nDimensions,nParticles,charge,alpha,beta);

#else

  Orbitals* dr = new Gaussians(nDimensions,nParticles,charge,alpha,beta);

#endif

#else

  Orbitals* dr = new  Diatomic(nDimensions,nParticles,charge,alpha,beta);

#endif

  dr->setR(R);

  wf = new WaveFunction(nParticles,nDimensions,alpha,beta,dr);

  LocalEnergyOpt* leo = new LocalEnergyOpt(nDimensions,nParticles,R,charge,alpha,beta,dr,wf);

  LocalEnergy* le = new LocalEnergy(nDimensions,nParticles,charge,alpha,beta,wf);
  LocalEnergyNumeric* len = new LocalEnergyNumeric(nDimensions,nParticles,charge,alpha,beta,wf);

  QuantumForce* qf = new QuantumForce(nDimensions,nParticles,charge,alpha,beta,dr,wf);

  QuantumForceNumeric* qfn = new QuantumForceNumeric(nDimensions,nParticles,charge,alpha,beta,wf);


  int accetpted = 0;

  double energySum = 0;
  double energySquaredSum = 0;

  double deltaE;

  double Rs;

  rOld = zeros( nParticles , nDimensions );
  rNew = zeros( nParticles , nDimensions );
 
  mat qforce_old  = zeros( nParticles , nDimensions );
  mat qforce_new  = zeros( nParticles , nDimensions );

  double ratio = 0;

  while (ratio == 0) {

    for ( int i = 0; i < nParticles ; i ++) {
      for ( int j = 0; j < nDimensions ; j ++) {
	rOld (i , j) = stepLength * ( ran2 (&idum) - 0.5) ;
      }
    } 

    rNew = rOld;

    wf->set_inverse(rOld);

    Rs = wf->slaterOpt(rNew,rOld,0);

    ratio = Rs*Rs*wf->jastrowOpt(rNew,rOld,0);
 
  }


  wf->set_inverse(rOld);

  //qfn->quantumforce(rOld , qforce_old,waveFunctionOld);

  qf->quantumforceOpt(rOld , qforce_old);


  for(int cycle = 0; cycle < nCycles; cycle++) {

    // New position to test
    for(int i = 0; i < nParticles; i++) {
      
      for(int j = 0; j < nDimensions; j++) {
	//rNew(i,j) = rOld(i,j) + stepLength*(ran2(&idum) - 0.5);
      	rNew(i,j) = rOld(i,j) + randn() * sqrt(timestep) + qforce_old(i, j) * timestep *D;
      }
      
      Rs = wf->slaterOpt(rNew,rOld,i);
      
      mat old_inv_up = wf->inv_up;
      mat old_inv_down = wf->inv_down;

      wf->update_inverse(Rs,rNew,i);

      //qfn->quantumforce( rNew , qforce_new,waveFunctionNew);
      qf->quantumforceOpt( rNew , qforce_new);
      
      double greensfunction = getGreensFunctionRatio(rNew,rOld,qforce_new,qforce_old);
      
      //       // Metropolis test
      if (ran2(&idum) <= greensfunction*Rs*Rs*wf->jastrowOpt(rNew,rOld,i)) {
      //if (ran2(&idum) <= wf->slaterbrute(rNew,rOld)*wf->jastrowbrute(rNew,rOld)) {
	for(int j = 0; j < nDimensions; j++) {
	  rOld(i,j) = rNew(i,j);
	  
	}
	
	qforce_old = qforce_new;

	accetpted++;
	
      } else {
	
	for(int j = 0; j < nDimensions; j++) {
	  rNew(i,j) = rOld(i,j);
	}

	wf->inv_up = old_inv_up;
	wf->inv_down = old_inv_down;

      }
      

      // update energies
      // deltaE = len->get_local_energy(rNew);
      // deltaE = le->get_local_energy(rNew);
      
#ifndef DIATOMIC

      deltaE = leo->get_local_energy(rNew);

#else

      deltaE = leo->get_local_energy_diatomic(rNew);

#endif

      energySum += deltaE;
      energySquaredSum += deltaE*deltaE;
    }
  }
  
  double energy = energySum/(nCycles * nParticles);
  double energySquared = energySquaredSum/(nCycles * nParticles);
  cout << "Energy: " << energy << " Energy (squared sum): " << energySquared << endl;
  cout << "Acceptred steps: " << (double)accetpted/(double)(nCycles * nParticles) << endl;        
  
  return energy;

}


double MCIntegrator::getGreensFunctionRatio(const mat &y, const mat &x, const mat &quantumForceNew, const mat &quantumForceOld){
  double argument1 = 0;
  double argument2 = 0;
  double argumentSum = 0;
  double greensFunctionRatio = 0;
  for (int i = 0; i < nParticles; i++){
    for (int j = 0; j < nDimensions; j++){
      argument1 += (y(i,j) - x(i,j) - 0.5*timestep*quantumForceOld(i,j))*(y(i,j) - x(i,j) - 0.5*timestep*quantumForceOld(i,j));
      argument2 += (x(i,j) - y(i,j) - 0.5*timestep*quantumForceNew(i,j))*(x(i,j) - y(i,j) - 0.5*timestep*quantumForceNew(i,j));
    }
  }
  argumentSum = (argument1 - argument2)/(2*timestep);
  greensFunctionRatio = exp(argumentSum);
  return greensFunctionRatio;
}

