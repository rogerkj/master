#include "mcintegrator.h"
#include "Minimizer.h"


#include <iostream>

#include <mpi.h>


using namespace std;

int main (int argc, char* argv[])
{
 
  int itermax = 100;

  int my_rank,num_procs;
 
  MPI_Status status;

  //Mpi initialization
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &num_procs);

 
  ofstream myfile;

  if (my_rank == 0) 
    myfile.open ("plot.dat");

  double left,rigth;
  
  double beta = 2;

  double h = 0.2;

  double k = pow(1.2,40);

  MCIntegrator *integrator = new MCIntegrator();

  for(int i = 0 ; i < itermax; i++) {

    double e = integrator->runMCIntegration(0,beta-h,0,my_rank);

    MPI_Reduce(&e,&left,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

    e = integrator->runMCIntegration(0,beta+h,0,my_rank);

    MPI_Reduce(&e,&rigth,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

    if (my_rank == 0) { 
      double delta = k*(rigth - left) /(2*pow(1.2,i+20)*num_procs);

      beta -= delta;
    
      myfile << (rigth+left)/(2*num_procs) << " " << beta << " " << delta << endl;

    }

     MPI_Bcast(&beta,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

  }

  if (my_rank == 0) 
    myfile.close();
    
  MPI_Finalize ();

  return(0);

}
