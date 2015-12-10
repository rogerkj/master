/*  Variational Quantum Montecarlo main algorithm
//
//
*/


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
  nDimensions(3),  // Number of dimensions
  charge(CHARGE), // Charge
  stepLength(0.1), //Steplength
  nParticles(NUMPART), //Number of particles
  idum(-1), //Random seed
  nCycles(NUMCYCLES),//Total number of cycles in MC
  timestep(0.005), // Timestep
    D(0.5),
  R(1.4) // Radius
   
{
}


double MCIntegrator::runMCIntegration(double a,double b,double _r,long _idum) {

  R = _r;

  alpha = a;
  beta = b;

  idum =  _idum;

  nThermalize = nCycles/10;

#ifndef DIATOMIC

#ifndef GAUSSIAN 

  //Orbital with hydrogenic basefunctions
  Orbitals* dr = new Hydrogenic(nDimensions,nParticles,charge,alpha,beta);

#else

  //Orbital with Gaussian basefunctions
  Orbitals* dr = new Gaussians(nDimensions,nParticles,charge,alpha,beta);

#endif

#else

  //Orbitals with diatomic basefunctions
  Orbitals* dr = new  Diatomic(nDimensions,nParticles,charge,alpha,beta);

#endif

  dr->setR(R);

  //Some tools 
  wf = new WaveFunction(nParticles,nDimensions,alpha,beta,dr);

  //Local energy code  (Optimized)
  LocalEnergyOpt* leo = new LocalEnergyOpt(nDimensions,nParticles,R,charge,alpha,beta,dr,wf);

  //Local energy code (Not optimized)
 // LocalEnergy* le = new LocalEnergy(nDimensions,nParticles,charge,alpha,beta,wf);

  //Local energy Numeric
  //LocalEnergyNumeric* len = new LocalEnergyNumeric(nDimensions,nParticles,charge,alpha,beta,wf);


  //Quantum force code
  QuantumForce* qf = new QuantumForce(nDimensions,nParticles,charge,alpha,beta,dr,wf);


  //Quantum force code (Numeric)
  //QuantumForceNumeric* qfn = new QuantumForceNumeric(nDimensions,nParticles,charge,alpha,beta,wf);

  //Keeps track of number of accepted steps
  int accetpted = 0;

  //Keeps track of energy and energy squared
  double energySum = 0;
  double energySquaredSum = 0;

  double deltaE;

  double Rs;

  //Positions of system
  rOld = zeros( nParticles , nDimensions );
  rNew = zeros( nParticles , nDimensions );
 
  //Quantum force matrises
  mat qforce_old  = zeros( nParticles , nDimensions );
  mat qforce_new  = zeros( nParticles , nDimensions );



  cout << "Starting vmc calculations." << endl;


  //This alogrithm picks out random starting positions
  //And checks if slaterdeterminant is too close to eachother
  mat rCurr;

  long oldidum = -1;

  bool not_done = true;
  while (not_done) {
    
    not_done = false;
    
    for(int i = 0; i < nParticles; i++) {
      for(int j = 0; j < nDimensions; j++) {
         rOld(i,j) = gaussianDeviate(&idum)*sqrt(timestep);
       }
    }
    
    rCurr = rOld;
    rNew = rOld;

    oldidum = idum;


    wf->set_inverse(rOld);
    
    qforce_old = qf->quantumforceOpt(rOld);

    for (int i = 0; i < nParticles; i++) {

      for(int j = 0; j < nDimensions; j++) {
         rNew(i,j) = rOld(i,j) + gaussianDeviate(&idum) * sqrt(timestep) + qforce_old(i, j) * timestep *D;
      }
      


      Rs = wf->slaterOpt(rNew,rOld,i);
      

      wf->update_inverse(Rs,rNew,i);

      cout << "Rs start: " << Rs << endl;

      if ((Rs == 0) || (isnan(Rs))) {
        not_done = true;

        cout << "slater wrong" << endl;

      }
      
      qforce_old = qf->quantumforceOpt(rNew);


      //rOld = rNew;

      rNew = rOld;

    }


    mat iup = wf->slater_up(rCurr);
    mat idown = wf->slater_down(rCurr);

    //      cout << inv(iup) << endl;
    //cout << inv(idown) << endl;

    for (int i = 0; i < nParticles/2; i++) {
        for (int j = 0; j < nParticles/2; j++) {

            if (i != j) {

                //Find distance between points
                double val_up = norm(iup.row(i) - iup.row(j));
                double val_down = norm(idown.row(i) - idown.row(j));

                // cout << val_up << "   "  << val_down << endl;

                //Check if distanses are too small
                if ((val_up < 0.1) || (val_down < 0.1)) {

                    not_done = true;

                }
            }
        }
    }



    for (int j = 0; j < nParticles; j++) {

      double val = norm(rCurr.row(j));

      if ((val < 0.1*sqrt(timestep)) || (val > 10*sqrt(timestep))) {

        not_done = true;

        cout << "center error" << endl;

      }



    }



    for (int i = 0; i < nParticles; i++) {
        for (int j = 0; j < nParticles; j++) {

            if (i != j) {

                 double val = norm(rCurr.row(i) - rCurr.row(j));

                  if (val < 0.1*sqrt(timestep)) {

                       not_done = true;

                       cout << "Electons too close" << endl;
                  }


            }
        }
    }

  }


  rOld = rCurr;
  rNew = rCurr;



idum = oldidum;


  wf->set_inverse(rOld);
  qforce_old = qf->quantumforceOpt(rOld);


  //Loop trough all cycles
  for(int cycle = 0; cycle < nCycles+nThermalize; cycle++) {



    // New position to test
    for(int i = 0; i < nParticles; i++) {
      

      //Make next move
      for(int j = 0; j < nDimensions; j++) {

	
            rNew(i,j) = rOld(i,j) + gaussianDeviate(&idum) * sqrt(timestep) + qforce_old(i, j) * timestep *D;
	
      }
      
      //Find slater determinant
      Rs = wf->slaterOpt(rNew,rOld,i);


      if ((Rs == 0) || (isnan(Rs))) {
        cout << "Rs = 0 or nan! in cycle " << cycle << " Particle nr: " << i << endl;

      }


      //Store current slaterdeterminant
      old_inv_up = wf->inv_up;
      old_inv_down = wf->inv_down;

      wf->update_inverse(Rs,rNew,i);
 
      //Find inverse
      //wf->set_inverse(rNew);


      //Find quantum force
      qforce_new = qf->quantumforceOpt( rNew );
      
      //get Greensfunction
      double greensfunction = getGreensFunctionRatio(rNew,rOld,qforce_new,qforce_old);

/*
      double greensfunction = 0.0;
      for (int k = 0; k < nParticles; k++)
      {
        for (int j = 0; j < nDimensions; j++)
        {
             greensfunction += (qforce_old(k,j) + qforce_new(k,j))*
                  (0.5*timestep*(qforce_old(k,j) - qforce_new(k,j)) + rOld(k,j) - rNew(k,j));
         }
       }

      greensfunction = exp(0.5*greensfunction);
*/
#ifdef JASTROW

      double jastrow = wf->jastrowOpt(rNew,rOld,i);

#else

      double jastrow = 1.0;

#endif

      
      // Metropolis test
      if (ran2(&idum) <= greensfunction*Rs*Rs*jastrow) {
      //if (ran2(&idum) <= wf->slaterbrute(rNew,rOld)*wf->jastrowbrute(rNew,rOld)) {

	//Step accepted,..   store result

	for(int j = 0; j < nDimensions; j++) {
	  rOld(i,j) = rNew(i,j);
	  
	}
	
	qforce_old = qforce_new;

	if (cycle > nThermalize)
	  accetpted++;
	
	//cout << "new step" << endl;

      } else {
	
	//Step not accepted... 

	for(int j = 0; j < nDimensions; j++) {
	  rNew(i,j) = rOld(i,j);
	}
	qforce_new = qforce_old;

	wf->inv_up = old_inv_up;
	wf->inv_down = old_inv_down;

    }

    if (cycle >= nThermalize) {
 
#ifndef DIATOMIC

	 //gets local energy for single atom 
	 deltaE = leo->get_local_energy(rNew);

#else

      //Gets local energy for diatiomic molecule
      deltaE = leo->get_local_energy_diatomic(rNew);

#endif

       //Sum energy and energy squared
       energySum += deltaE;
       energySquaredSum += deltaE*deltaE;
    

     }
   }
  }

  //Get average energy and energy squared
  double energy = energySum/(nCycles * nParticles);
  double energySquared = energySquaredSum/(nCycles * nParticles);


  //Print result to screen
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

