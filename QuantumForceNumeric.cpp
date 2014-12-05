#include "WaveFunction.h"
#include "QuantumForceNumeric.h"
#include "lib.h"

#include <armadillo>
#include <iostream>

using namespace arma;
using namespace std;

QuantumForceNumeric::QuantumForceNumeric(int nDim,int nPart,int ch,double a,double b,WaveFunction* _wf) :
  h(0.001),
  h2(1000000)

{
	nDimensions = nDim;
	nParticles = nPart;
	charge = ch;
	alpha = a;
	beta = b;

	wf = _wf;
}

void QuantumForceNumeric::quantumforce (const mat &r , mat &qforce ,double wfunc ) {

  double wfminus , wfplus ;
  mat r_plus  = zeros<mat>( nParticles , nDimensions );
  mat r_minus = zeros<mat>( nParticles , nDimensions );

  r_plus = r_minus = r;

  // compute the firstderivative
  for ( int i = 0; i <  nParticles ; i ++) {
    for ( int j = 0; j <  nDimensions ; j ++) {
      r_plus(i, j) =  r(i, j) + h;
      r_minus(i, j) = r(i, j) - h;

      wfminus =  wf->slater(r_minus);
      wfplus =  wf->slater(r_plus);

      qforce(i, j) = (wfplus - wfminus) * 2.0 / wfunc / (2.0 * h);
    
      r_plus(i, j) = r(i, j) ;
      r_minus( i, j) = r(i, j) ;
    }
  }
} // end of quantum force function

