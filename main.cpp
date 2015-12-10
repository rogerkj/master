#include "mcintegrator.h"
#include "Minimizer.h"


#include <iostream>

#include <mpi.h>


using namespace std;

int main (int argc, char* argv[])
{
 

  int my_rank,num_procs;
 
  MPI_Status status;

  //Mpi initialization
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &num_procs);

 
  int totn = 1000;

  int n = 1;//totn/num_procs;
  

  ofstream myfile;

  double* res;


  if (my_rank == 0) {
    myfile.open ("plot.dat");

    res = (double*) malloc(n*num_procs*3*sizeof(double));

  }

  double* e = (double*) malloc(3*sizeof(double));

  
  MCIntegrator *integrator = new MCIntegrator();

  for (int i = 0; i < n; i++) {

    double alpha = 0;
    
    double beta = 1.66201;//0.1 + 7.0*(num_procs*i + my_rank)/(n*num_procs);
    //   double beta = 0;//0.1 + 3.0*(num_procs*i + my_rank)/(n*num_procs);
    
    e[0] = integrator->runMCIntegration(alpha,beta,0,-3011);//-my_rank-i*num_procs);
    e[1] = beta;
    e[2] = -my_rank - i*num_procs;
    
    
    MPI_Gather(e,3,MPI_DOUBLE,&res[i*num_procs*3],3,MPI_DOUBLE,0,MPI_COMM_WORLD);
    
    
  }
  
  if (my_rank == 0) { 
    
    for (int t = 0 ; t < n*num_procs; t++) {
	
      myfile << res[t*3] << " " << res[t*3+1] << " " << res[t*3+2] << endl;
      
    }
    
    myfile.close();
     
  }
  
  MPI_Finalize();
  
  return(0);
}
