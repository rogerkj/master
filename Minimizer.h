#ifndef MINIMIZER_H
#define MINIMIZER_H

#include <armadillo>
#include "mcintegrator.h"

using namespace arma;


class Minimizer
{

 public:
  
  Minimizer();

  double find_minimum_beta(MCIntegrator* mci,double _beta,double _r);
  double find_minimum_beta2(MCIntegrator* mci,double _beta,double _r);


 private:


};

#endif
