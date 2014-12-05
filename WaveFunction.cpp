#include <armadillo>
#include <iostream>

#include "WaveFunction.h"

using namespace arma;
using namespace std;

WaveFunction::WaveFunction(int nPart,int nDim,double a,double b,Orbitals* _dr) {

  nDimensions = nDim;
  nParticles = nPart;

  alpha = a;
  beta = b;

  //  amat  = zeros<mat>( 4 , 4 );

  spin = new bool[nParticles];

  for (int i = 0; i < nParticles; i++) {
    if (i < nParticles/2)
      spin[i] = true;
    else
      spin[i] = false;
  }
  
  dr = _dr;
 
}


void WaveFunction::update_inverse(double Rs,const mat &r,int i) {
  
  vec I(nParticles);

  if(i < nParticles/2) {

    for(int j = 0; j < nParticles/2; j++) {
      I(j) = 0;
      if(i!=j) {
	for(int l = 0; l < nParticles/2; l++)
	  I(j) +=  dr->waveFunction(r,i,l)*inv_up(l,j);
      }
    }
    for(int j = 0; j < nParticles/2; j++) {
      if (i!=j) {
	for(int k = 0; k < nParticles/2; k++)
	  inv_up(k,j) = inv_up(k,j) - I(j)*inv_up(k,i)/Rs;
      }
    }

    for (int k = 0; k < nParticles/2; k++)
      inv_up(k,i) = inv_up(k,i)/Rs;

  } else {

    i = i - nParticles/2;

   for(int j = 0; j < nParticles/2; j++) {
      I(j) = 0;
      if(i!=j) {
	for(int l = 0; l < nParticles/2; l++)
	  I(j) +=  dr->waveFunction(r,i+nParticles/2,l)*inv_down(l,j);
      }
    }
    for(int j = 0; j < nParticles/2; j++) {
      if (i!=j) {
	for(int k = 0; k < nParticles/2; k++)
	  inv_down(k,j) = inv_down(k,j) - I(j)*inv_down(k,i)/Rs;
      }
    }

    for (int k = 0; k < nParticles/2; k++)
      inv_down(k,i) = inv_down(k,i)/Rs;

  }
  
}

void WaveFunction::set_inverse(const mat &r) {

  inv_up = inv(slater_up(r));
  inv_down = inv(slater_down(r));

}

double WaveFunction::slaterOpt(const mat &rNew,const mat &rOld,int i) {
    
  double Rs = 0;
  
  if (i < nParticles/2) {
    for (int j = 0 ; j < nParticles/2; j++) {
      Rs += dr->waveFunction(rNew,i,j)*inv_up(j,i);
    }
  } else {
    for (int j = 0 ; j < nParticles/2; j++) {
      Rs += dr->waveFunction(rNew,i,j)*inv_down(j,i-nParticles/2);
    }
  }
  
  return Rs;
}


double WaveFunction::jastrowOpt(const mat &rNew, const mat &rOld,int i) {

  double sum = 0;

  for (int j = 0; j < nParticles; j++) {

    rowvec r12_new = rNew.row(i) - rNew.row(j);
    double r12norm_new = norm(r12_new, 2);
  
    rowvec r12_old = rOld.row(i) - rOld.row(j);
    double r12norm_old = norm(r12_old, 2);

    sum += amat(i,j)*(r12norm_new / (1 + beta * r12norm_new) - r12norm_old / (1 + beta * r12norm_old));

  }

  double e = exp(sum);

  return e*e;

}

double WaveFunction::amat(int i, int j) {
  
  if (spin[i] == spin[j])
    return 0.25;
  else
    return 0.5;
  
}
     


double WaveFunction::jastrowbrute(const mat &r_new,const mat &r_old) {

  double jast_new = 1.0;
  double jast_old = 1.0;
  
  for(int i = 0; i < nParticles; i++) {
    for(int j = i+1; j < nParticles; j++) {
      jast_new *= jastrowFactor(r_new,i,j);
      jast_old *= jastrowFactor(r_old,i,j);
    }
  }

  return jast_new*jast_new/(jast_old*jast_old);

}


mat WaveFunction::slater_up(const mat &r) {
  
  mat up_mat  = zeros<mat>( nParticles/2 , nParticles/2 );

  for(int i = 0;i < nParticles/2; i++) {
    for(int j = 0; j < nParticles/2; j++) {
      up_mat(i,j)  = dr->waveFunction(r,i,j);
    }
  }

  return up_mat;

}

mat WaveFunction::slater_down(const mat &r) {
  
  mat down_mat  = zeros<mat>( nParticles/2 , nParticles/2 );

  for(int i = 0;i < nParticles/2; i++) {
    for(int j = 0; j < nParticles/2; j++) {
      down_mat(i,j)  = dr->waveFunction(r,i+nParticles/2,j);
    }
  }

  return down_mat;

}


double WaveFunction::slaterbrute(const mat &rNew,const mat &rOld) {

  mat det_up_new = slater_up(rNew);
  mat det_down_new = slater_down(rNew);

  mat det_up_old = slater_up(rOld);
  mat det_down_old = slater_down(rOld);

  double detnew = det(det_up_new)*det(det_down_new);
  double detold = det(det_up_old)*det(det_down_old);
  

  return detnew*detnew/(detold*detold);
}

double WaveFunction::slater(const mat &r) {

  mat detPlus  = zeros<mat>( nParticles/2 , nParticles/2 );
  mat detMinus = zeros<mat>( nParticles/2 , nParticles/2 );

  for(int i = 0;i < nParticles/2; i++) {
    for(int j = 0; j < nParticles/2; j++) {
      detPlus(i,j)  = dr->waveFunction(r,i,j);
      detMinus(i,j) = dr->waveFunction(r,i+nParticles/2,j);
    }
  }

  double jast = 1.0;
  
  for(int i = 0; i < nParticles; i++) {
    for(int j = i+1; j < nParticles; j++) {
      jast *= jastrowFactor(r,i,j);
    }
  }

  return det(detPlus)*det(detMinus) * jast;

}
/*

double WaveFunction::waveFunction(const mat &r,int nParticle,int orbital) {
  
  int i = nParticle;
 
  double argument = norm(r.row(i));
 
  switch (orbital) {
  case 0:
    return exp(-argument * alpha);
  case 1:
    return (alpha*argument - 2)*exp(-argument * alpha/2.0);
  case 2: 
    return r(i,0) * exp(-alpha * argument / 2.0); 
  case 3: 
    return r(i,1) * exp(-alpha * argument / 2.0); 
  case 4: 
    return r(i,2) * exp(-alpha * argument / 2.0); 

  }

  return 0;

}
*/

double WaveFunction::jastrowFactor(const mat &r,int i,int j) {

  rowvec r12 = r.row(i) - r.row(j);
  double r12norm = norm(r12);
  
  return exp(amat(i,j)*r12norm / (1 + beta * r12norm));
}
