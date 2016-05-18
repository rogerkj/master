#include "Minimizer.h"



#include <armadillo>
#include <iostream>

using namespace arma;
using namespace std;

Minimizer::Minimizer() {}

double Minimizer::find_minimum_beta2(MCIntegrator* mci,double _beta,double _r) {

  int itermax = 200;

  double step_length = 0.35;
  double beta = _beta;

  double left,rigth;


  int m = 10;

  //  ofstream myfile;
  //myfile.open ("plot.dat");

  for(int i = 50 ; i < itermax; i++) {

    cout << "iteration numner: " << i << endl;

    double gradient = 0;

    for (int j = 0; j < m; j++) {
   
      rigth =  mci->runMCIntegration(0,beta-step_length,_r,-j-i-1);
      left =  mci->runMCIntegration(0,beta+step_length,_r,-j-i-1);

      gradient += left - rigth;

    }

    double scaling = (double) i;

    double forskjell = (double) gradient / pow(scaling,1.2);

    forskjell *= 50;
      
    beta -= forskjell;

    cout << "Beta: " << beta << endl;

    //myfile << i << " " << beta << " " << rigth << endl;

  }

  // myfile.close();

  return beta;


}

double Minimizer::find_minimum_beta(MCIntegrator* mci,double _beta,double _r) {

  int itermax = 200;

  double step_length = 0.2;
  double beta = _beta;

  int m = 10;

  //ofstream myfile;
  //myfile.open ("plot.dat");


  double sga_upp;
  double sga_org;

  double sga_upp_sum;
  double sga_org_sum;

  double sga_tot_sum;

  double w_i;

  double m_tilde;

  double abvar;

  for(int i = 50 ; i < itermax; i++) {

    cout << "iteration numner: " << i << endl;

    m_tilde = 0;
    sga_upp = 0;
    sga_org = 0;

    sga_upp_sum = 0;
    sga_org_sum = 0;
    

    for (int j = 0; j < m; j++) {
   
      sga_org =  mci->runMCIntegration(0,beta-step_length,_r,-j-i-1);
      sga_upp =  mci->runMCIntegration(0,beta+step_length,_r,-j-i-1);

      w_i = sga_upp*sga_upp / ( sga_org * sga_org);

      m_tilde += w_i;

      sga_org_sum += w_i * sga_org;

      sga_upp_sum += w_i *  sga_upp;
    }

    sga_tot_sum = sga_upp_sum - sga_org_sum;
    sga_tot_sum = sga_tot_sum / (m_tilde*step_length*2);

    double scaling = (double) i;

    double forskjell = (double) sga_tot_sum / scaling;

    forskjell *= 500;
      

    beta -= forskjell;

    cout << "Beta: " << beta << endl;

    // myfile << i << " " << beta << endl;

  }

  //myfile.close();
 
  return beta;
 

}
