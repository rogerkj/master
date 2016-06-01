#include "mcintegrator.h"
#include "setup.h"

#include <iostream>

using namespace std;

struct res {

  double E;
  int blockSize;
  double sigma;

};


int main (int argc, char* argv[])
{

  double alpha = ALPHA;
  double beta = BETA;

  double R = RADIUS;



  int deltaBlockSize = 10;

  int maxBlockSizeTreshold = 1E4;

  int stepLength = 10;


 

  MCIntegrator *integrator = new MCIntegrator();

  retvalues rv = integrator->runMCIntegration(alpha,beta,R,-1);
    
  int N = integrator->nCycles;

  int maxBlockSize = (int) N / 10;


  if(maxBlockSize > maxBlockSizeTreshold)
    maxBlockSize = maxBlockSizeTreshold;

  res* block = (res*) malloc(sizeof(res)*maxBlockSize);///deltaBlockSize);

  cout << rv.energy << " " << rv.energySquared << endl;

  cout << ((rv.energy*rv.energy) - rv.energySquared)  << endl;

  cout << "Doing blocking." << endl;

  for (int n = 1 ; n*deltaBlockSize <= maxBlockSize; n++) {

    int blockSize = n*deltaBlockSize;

    int samples = integrator->nCycles;

    int blocks = int(samples/blockSize);


    double* averageEnergyBlock = (double*) malloc(sizeof(double)*blocks);

    memset(averageEnergyBlock,0,sizeof(double)*blocks);


    double energySum;


    for (int i = 0; i < blocks; i++) {
	
      energySum = 0;


        for (int j = i * blockSize; j < i * blockSize + blockSize; j++){
            energySum += rv.local_energy[j];
        }
        averageEnergyBlock[i] = energySum / blockSize;
    }

    // Calculating the mean and variance of all the blocks
    double E = 0;
    double E2 = 0;

    for (int i = 0; i < blocks; i++) {
        E += averageEnergyBlock[i];
        E2 += averageEnergyBlock[i] * averageEnergyBlock[i];
    }

    // Averaging
    E  /= blocks;
    E2 /= blocks;

    double sigma = E2 - E*E;
    sigma = sqrt(sigma/blocks);

    block[n-1].blockSize = blockSize;
    block[n-1].E = E;
    block[n-1].sigma = sigma;

    free(averageEnergyBlock);

  }

  
  cout << "Writing to file." << endl;

  ofstream myfile;
  myfile.open ("blocks.dat");

  myfile << "#VMC blocks\n";
 
  for (int t = 1; t < maxBlockSize; t++) {

    myfile << block[t].blockSize << "\t" << block[t].sigma << "\t" << block[t].E << "\n";

  }
    

  myfile.close();
   
  return(0);

}
