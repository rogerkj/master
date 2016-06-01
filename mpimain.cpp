#include "setup.h"
#include "mcintegrator.h"
#include "Minimizer.h"


#include <iostream>

#include <mpi.h>


using namespace std;

int main (int argc, char* argv[])
{
 

  double R = 2.282;

  double alpha = CHARGE;
  double beta = 0.3;


  double delta_alpha = 0.1;
  double delta_beta = 0.03;


  int my_rank,num_procs;
 
  MPI_Status status;

  //Mpi initialization
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &num_procs);



  int n = 4000;//totn/num_procs;
  

  double* res;

  ofstream myfile;


  if (my_rank == 0) {

    myfile.open ("plot.dat");

    res = (double*) malloc(n*3*sizeof(double));

  }

  double e;

  double left_alpha,right_alpha;
  double left_beta,right_beta;
  double middle;


  MCIntegrator *integrator = new MCIntegrator();

  for (int i = 0; i < n; i++) {
     
    e = integrator->runMCIntegration(alpha,beta,R,-(i*num_procs + my_rank)*5);
    MPI_Reduce(&e,&middle,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

    e = integrator->runMCIntegration(alpha-delta_alpha,beta,R,-(i*num_procs + my_rank)*5-1);
    MPI_Reduce(&e,&left_alpha,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    
    e = integrator->runMCIntegration(alpha+delta_alpha,beta,R,-(i*num_procs + my_rank)*5-2);
    MPI_Reduce(&e,&right_alpha,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

    e = integrator->runMCIntegration(alpha,beta-delta_beta,R,-(i*num_procs + my_rank)*5-3);
    MPI_Reduce(&e,&left_beta,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    
    e = integrator->runMCIntegration(alpha,beta+delta_beta,R,-(i*num_procs + my_rank)*5-4);
    MPI_Reduce(&e,&right_beta,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);



    if (my_rank == 0) {
    
      res[i*3] = middle / num_procs;
      res[i*3+1] = alpha;
      res[i*3+2] = beta;

      left_alpha = left_alpha/num_procs;
      right_alpha = right_alpha/num_procs;

      left_beta = left_beta/num_procs;
      right_beta = right_beta/num_procs;

      middle = middle / num_procs;


      double dalpha = (right_alpha - left_alpha) / (2*delta_alpha);
      double dbeta  = (right_beta - left_beta) / (2*delta_beta);

      double d2alpha = (right_alpha - 2*middle + left_alpha)/(delta_alpha*delta_alpha);
      double d2beta = (right_beta - 2*middle + left_beta)/(delta_beta*delta_beta);

      alpha -= dalpha / d2alpha;


      double newbeta = beta - dbeta / d2beta;

      if (newbeta > 0)
	beta = newbeta;


      myfile << res[i*3] << " " << res[i*3+1] << " " << res[i*3+2] << " " << i*5 << endl;

    }

    MPI_Bcast(&alpha,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&beta,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    
      
  }
   
  if (my_rank == 0) { 
    /*  
    for (int t = 0 ; t < n; t++) {
	
     
      
    }
    */
    myfile.close();
     
  }
  
  MPI_Finalize();
  
  return(0);
}
