#ifndef LOSSFUNC_H_
#define LOSSFUNC_H_

#include<cmath>

inline double d_L_mu(double y, double mu){
  return -y * 1.0/(1.0+exp(y * mu));
}

inline double d2_L_mu(double y, double mu){
  return (y * y) * 1.0/(1.0+exp(-y * mu)) * 1.0/(1.0+exp(y * mu));
  //double pexp = exp(-y * mu);
  //double plog = 1.0+pexp;
  //return (y * y) *  pexp / (plog*plog);
}

inline double L(double y, double mu){
  return log(1.0+exp(-y*mu));
}

#endif /*LOSSFUNC_H_*/
