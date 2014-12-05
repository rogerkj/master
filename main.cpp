#include "mcintegrator.h"

#include <iostream>

using namespace std;

int main (int argc, char* argv[])
{

  MCIntegrator *integrator = new MCIntegrator();
   
  
  int N = 100;
  double db = (1.0)/N;

  double b = 0.01;
  
 
  
  
  ofstream myfile;
  myfile.open ("plot.dat");

  myfile << "#VMC plot\n";
 
  for (int t = 0; t < N ; t++) {
  
    double e = integrator->runMCIntegration(0,b);
    
    myfile << b << "\t" << e << "\n";
    b += db;

  }
    

  myfile.close();
    
    
  return(0);

}
