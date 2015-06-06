#include "pulseShape.hh"
#include <math.h>

double pulseShape::eval(double x, float tZero, float alpha, float beta){
     float t = x - tZero;
     double tmp = 1. + t / (alpha * beta);
     if(tmp<1e-9){
	  return 0.;
     }else{
	  return pow( tmp, alpha ) * exp( -t / beta );
     }
}

#ifdef shervin
void pulseShape::NoiseCorrInit()
{
  // covariance
  double Ctmp[MAX_SAMPLES];
  for(int i=0; i<MAX_SAMPLES; ++i) Ctmp[i] = EECorrNoiseMatrixG12[i];
  
  // reset
  double L[MAX_SAMPLES][MAX_SAMPLES];
  for(int i=0; i<MAX_SAMPLES; ++i){
    for(int j=0; j<MAX_SAMPLES; ++j){
      L[i][j]=0;
    }
  }
  // decomposition
  L[0][0] = sqrt(Ctmp[0]);
  for( int col=1; col<MAX_SAMPLES; col++){
    L[0][col]=0;
  }
  for( int row=1; row<MAX_SAMPLES; row++){
    for( int col=0; col<row; col++ ){
      double sum1 = 0;
      int m=abs(row-col);
      for( int k=0; k<col; ++k) sum1 += L[row][k]*L[col][k];
      L[row][col] = (Ctmp[m] - sum1)/L[col][col];
    }
    double sum2 = 0;
    for( int k=0; k<row; ++k) sum2 += L[row][k]*L[row][k];
    L[row][row] = sqrt( Ctmp[0] - sum2 );
    for( int col=row+1; col<MAX_SAMPLES; col++ ) L[row][col] = 0;
  }
  // save result
  for( int row=0; row<MAX_SAMPLES; row++){
    for( int col=0; col<MAX_SAMPLES; col++ ){
      LGen[row][col] = L[row][col];
    }
  }
}
#endif
