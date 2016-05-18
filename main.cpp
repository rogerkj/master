#include "setup.h"
#include "mcintegrator.h"
#include "Minimizer.h"


#include <iostream>

#include <mpi.h>


using namespace std;

int main (int argc, char* argv[])
{

  double R = 2.3481;


  int my_rank,num_procs;
 
  MPI_Status status;

  //Mpi initialization
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &num_procs);



  ofstream myfile;

  if (my_rank == 0) {

    myfile.open ("plot.dat");

  }

  double e,sum;


  MCIntegrator *integrator = new MCIntegrator();

 
     
  e = integrator->runMCIntegration(0,0,R, -my_rank);
  MPI_Reduce(&e,&sum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);


  if (my_rank == 0) {
     
    myfile << sum/num_procs << endl;
    
    myfile.close();
     
  }
  
  MPI_Finalize();
  
  return(0);
}
