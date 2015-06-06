#include "TRandom3.h"

TRandom rnd;
const double tPhase = 0;

double funPulseShape(double *x, double *par)
{
  double t = x[0] - par[0];
  double alpha = par[1];
  double beta  = par[2];
  double tmp = 1. + t / (alpha * beta);
  if(tmp<1e-9){
    return 0.;
  }else{
    return pow( tmp, alpha ) * exp( -t / beta );
  }  
}
