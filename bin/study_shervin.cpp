//g++ -o study_Ra study_Ra.cc PulseChiSqSNNLS.cc -std=c++11 `root-config --cflags --glibs` 


#include <iostream>
#include "PulseChiSqSNNLS.h"
#include "MultiFitToyMC.h"

#include "Riostream.h"
#include "TTree.h"
#include "TF1.h"
#include "TProfile.h"
#include "TH2.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TFile.h"
#include "TStyle.h" 

#include "digi.hh"

#define NSAMPLES 10
#define NTOYS 1000
#define nOOT 14
#define startOOT -12
#define stopOOT 4
#define ADC2GeVEB 0.03
#define ADC2GeVEE 0.069

using namespace std;

double tZero = 0.0;
double alpha = 1.63669;
double beta  = 1.43412;

double weights[NSAMPLES] = {
  -0.409783, -0.409783, -0.409783, 0.0, 0.258171,
  0.318941, 0.283887, 0.216833, 0.151517, 0.0 
};


/*
double weights[NSAMPLES] = {
  -0.05, -0.15, -0.3, -0.625, -0.7, 0.375, 
  1., 1.43, 0., 0.};
*/

FullSampleVector fullpulse(FullSampleVector::Zero());
FullSampleMatrix fullpulsecov(FullSampleMatrix::Zero());
SampleMatrix noisecor(SampleMatrix::Zero());
BXVector activeBX;
SampleVector amplitudes(SampleVector::Zero());
SampleVector amplitudesNoise(SampleVector::Zero());
SampleVector amplitudesRef(SampleVector::Zero());

double EECorrNoiseMatrixG12[] = {
  1.00000, 0.71373, 0.44825, 0.30152, 0.21609,
  0.14786, 0.11772, 0.10165, 0.09465, 0.08098 
};
double LGen[NSAMPLES][NSAMPLES];

TF1 *fps;
TFile *fout;

/*
TProfile* hMF_ampRef_vs_ampTrue;
TH1F*     hMF_Delta_ampRef_ampNoise;
TProfile* hMF_relDelta_ampRef_amp_vs_BX;
TProfile* hMF_relDelta_ampNoise_amp_vs_BX;
TProfile* hMF_relDelta_ampTrue_amp_vs_BX;
TProfile* hMF_Ratio_ampRef_amp_vs_BX;
TProfile* hMF_Ratio_ampNoise_amp_vs_BX;
TProfile* hMF_Ratio_ampTrue_amp_vs_BX;



TProfile* hWg_ampRef_vs_ampTrue;
TH1F*     hWg_Delta_ampRef_ampNoise;
TProfile* hWg_relDelta_ampRef_amp_vs_BX;
TProfile* hWg_relDelta_ampNoise_amp_vs_BX;
TProfile* hWg_relDelta_ampTrue_amp_vs_BX;
TProfile* hWg_Ratio_ampRef_amp_vs_BX;
TProfile* hWg_Ratio_ampNoise_amp_vs_BX;
TProfile* hWg_Ratio_ampTrue_amp_vs_BX;
*/


TProfile* hMF_ampRef_vs_BX;
TProfile* hMF_ampNoise_vs_BX;
TProfile* hMF_ampPU_vs_BX;
TProfile* hWg_ampRef_vs_BX;
TProfile* hWg_ampNoise_vs_BX;
TProfile* hWg_ampPU_vs_BX;

void NoiseCorrInit()
{
  // covariance
  double Ctmp[NSAMPLES];
  for(int i=0; i<NSAMPLES; ++i) Ctmp[i] = EECorrNoiseMatrixG12[i];
  
  // reset
  double L[NSAMPLES][NSAMPLES];
  for(int i=0; i<NSAMPLES; ++i){
    for(int j=0; j<NSAMPLES; ++j){
      L[i][j]=0;
    }
  }
  // decomposition
  L[0][0] = sqrt(Ctmp[0]);
  for( int col=1; col<NSAMPLES; col++){
    L[0][col]=0;
  }
  for( int row=1; row<NSAMPLES; row++){
    for( int col=0; col<row; col++ ){
      double sum1 = 0;
      int m=abs(row-col);
      for( int k=0; k<col; ++k) sum1 += L[row][k]*L[col][k];
      L[row][col] = (Ctmp[m] - sum1)/L[col][col];
    }
    double sum2 = 0;
    for( int k=0; k<row; ++k) sum2 += L[row][k]*L[row][k];
    L[row][row] = sqrt( Ctmp[0] - sum2 );
    for( int col=row+1; col<NSAMPLES; col++ ) L[row][col] = 0;
  }
  // save result
  for( int row=0; row<NSAMPLES; row++){
    for( int col=0; col<NSAMPLES; col++ ){
      LGen[row][col] = L[row][col];
    }
  }
}


void initHist()
{
  gStyle->SetOptStat(1111);

  //  std::cout << "" << std::endl;
  fout = new TFile("output.root","recreate");

  /*
  hMF_ampRef_vs_ampTrue = new TProfile("hMF_ampRef_vs_ampTrue", "", 100, 0., 10.);
  hMF_Delta_ampRef_ampNoise = new TH1F("hMF_Delta_ampRef_ampNoise", "", 500, -5., 5.);

  hMF_relDelta_ampRef_amp_vs_BX = new TProfile("hMF_relDelta_ampRef_amp_vs_BX", "", 30, -15., 15.);
  hMF_relDelta_ampNoise_amp_vs_BX = new TProfile("hMF_relDelta_ampNoise_amp_vs_BX", "", 30, -15., 15.);
  hMF_relDelta_ampTrue_amp_vs_BX = new TProfile("hMF_reldelta_ampTrue_amp_vs_BX", "", 30, -15., 15.);
  hMF_Ratio_ampRef_amp_vs_BX = new TProfile("hMF_Ratio_ampRef_amp_vs_BX", "", 30, -15., 15.);
  hMF_Ratio_ampNoise_amp_vs_BX = new TProfile("hMF_Ratio_ampNoise_amp_vs_BX", "", 30, -15., 15.);
  hMF_Ratio_ampTrue_amp_vs_BX = new TProfile("hMF_Ratio_ampTrue_amp_vs_BX", "", 30, -15., 15.);


  hWg_ampRef_vs_ampTrue = new TProfile("hWg_ampRef_vs_ampTrue", "", 100, 0., 10.);
  hWg_Delta_ampRef_ampNoise = new TH1F("hWg_Delta_ampRef_ampNoise", "", 500, -5., 5.);

  hWg_relDelta_ampRef_amp_vs_BX = new TProfile("hWg_relDelta_ampRef_amp_vs_BX", "", 30, -15., 15.);
  hWg_relDelta_ampNoise_amp_vs_BX = new TProfile("hWg_relDelta_ampNoise_amp_vs_BX", "", 30, -15., 15.);
  hWg_relDelta_ampTrue_amp_vs_BX = new TProfile("hWg_reldelta_ampTrue_amp_vs_BX", "", 30, -15., 15.);
  hWg_Ratio_ampRef_amp_vs_BX = new TProfile("hWg_Ratio_ampRef_amp_vs_BX", "", 30, -15., 15.);
  hWg_Ratio_ampNoise_amp_vs_BX = new TProfile("hWg_Ratio_ampNoise_amp_vs_BX", "", 30, -15., 15.);
  hWg_Ratio_ampTrue_amp_vs_BX = new TProfile("hWg_Ratio_ampTrue_amp_vs_BX", "", 30, -15., 15.);
  */

  hMF_ampRef_vs_BX = new TProfile("hMF_ampRef_vs_BX", "", 30, -15., 15.);
  hMF_ampNoise_vs_BX = new TProfile("hMF_ampNoise_vs_BX", "", 30, -15., 15.);
  hMF_ampPU_vs_BX = new TProfile("hMF_ampPU_vs_BX", "", 30, -15., 15.);
  hWg_ampRef_vs_BX = new TProfile("hWg_ampRef_vs_BX", "", 30, -15., 15.);
  hWg_ampNoise_vs_BX = new TProfile("hWg_ampNoise_vs_BX", "", 30, -15., 15.);
  hWg_ampPU_vs_BX = new TProfile("hWg_ampPU_vs_BX", "", 30, -15., 15.);
  
  
}

void init()
{
  initHist();
  NoiseCorrInit();

  fps = new TF1("fps",funPulseShape, -10., 100., 3);
  fps->SetParameters(tZero+2.0, alpha, beta);
  double pulseShapeTemplate[12];
  for(int i=0; i<12; i++) pulseShapeTemplate[i] = fps->Eval((double)i);
  for(int i=0; i<12; i++) pulseShapeTemplate[i] /= pulseShapeTemplate[2];
  for (int i=0; i<12; ++i) fullpulse(i+7) = pulseShapeTemplate[i];

  for (int i=0; i<NSAMPLES; ++i) {
    for (int j=0; j<NSAMPLES; ++j) {
      int vidx = std::abs(j-i);
      noisecor(i,j) = EECorrNoiseMatrixG12[vidx];
    }
  }

  int activeBXs[] = { -5, -4, -3, -2, -1,  0,  1,  2,  3,  4 };
  activeBX.resize(NSAMPLES);
  for (unsigned int ibx=0; ibx<NSAMPLES; ++ibx) {
    activeBX.coeffRef(ibx) = activeBXs[ibx];
  } 
}


vector<int> newHit(double aMax, double trueRMS, double pedestal)
{
  double amp[NSAMPLES];
  for(int i=0; i<NSAMPLES; i++) amp[i] = 0.;

  fps->SetParameters(tZero + 5.0, alpha, beta);
  for(int i=0; i<NSAMPLES; i++) amp[i] += aMax * fps->Eval((double)i);

  for(int i=0; i<NSAMPLES; i++){
    amp[i] += pedestal;
  }

  // add noise
  double noiseA[NSAMPLES];
  for(int i=0; i<NSAMPLES; i++){
    noiseA[i] = rnd.Gaus(0.0, 1.0);
  }
  // correlations
  double noiseB[NSAMPLES];
  for(int i=0; i<NSAMPLES; ++i){
    noiseB[i]=0;
    for(int j=0; j<NSAMPLES; ++j) noiseB[i] += LGen[i][j]*noiseA[j];
  }
  for(int i=0; i<NSAMPLES; ++i){
    noiseA[i] = noiseB[i];
    amp[i] += trueRMS * noiseA[i];
  }

  vector<int> a;
  for(int i=0; i<NSAMPLES; i++) a.push_back((int)amp[i]);
  return a;
}
 

//vector<int> newDigi( const double *e, double trueRMS, double pedestal, double aMax, int tPU)
// this function builds the pulse shape for a PU contribution with energy e and bunch corss tPU with respect to the in-time bunch cross:
// 0 = in-time
// -1 = early OOT PU
// +1 = late OOT PU
vector<float> newDigiPU( double e, int tPU)
{
  //  std::cout << " >>> chiamata newDigiPU " << std::endl;
  //  building pulse shape to be fitted
  vector<float> ampPU;

  //add PU
  fps->SetParameters(tPU + 5.0, alpha, beta);
  for(int i=0; i<NSAMPLES; i++) {
    ampPU.push_back( e * fps->Eval((double)i) );
    //    std::cout << "ampPU["<<i<<"] = " << ampPU.at(i) << std::endl; 
  }
  
  return ampPU;

}


//vector<int> newDigi( const double *e, double trueRMS, double pedestal, double aMax, int tPU)
vector<float> newDigi( double e, double trueRMS, double pedestal, double aMax, vector<float> aMaxPU, int tPU)
{

  //     std::cout << " >>> chiamata  " << std::endl;
  //  building pulse shape to be fitted
  double amp[NSAMPLES];

  //add PU
  //   fps->SetParameters(tPU + 5.0, alpha, beta);
  for(int i=0; i<NSAMPLES; i++) amp[i] = aMaxPU.at(i);

  //add signal
  fps->SetParameters(tZero + 5.0, alpha, beta);
  for(int i=0; i<NSAMPLES; i++) {
    amp[i] += aMax * fps->Eval((double)i);
    //    std::cout << "amp["<<i<<"] = " << amp[i] << std::endl; 
  }

  //add pedestal
  if(pedestal > 0.){
    for(int i=0; i<NSAMPLES; i++){
      amp[i] += pedestal;
    }
  }

  // add noise
  if(trueRMS > 0.){
    double noiseA[NSAMPLES];
    for(int i=0; i<NSAMPLES; i++){
      noiseA[i] = rnd.Gaus(0.0, 1.0);
    }
    // correlations
    double noiseB[NSAMPLES];
    for(int i=0; i<NSAMPLES; ++i){
      noiseB[i]=0;
      for(int j=0; j<NSAMPLES; ++j) noiseB[i] += LGen[i][j]*noiseA[j];
    }
    for(int i=0; i<NSAMPLES; ++i){
      noiseA[i] = noiseB[i];
      amp[i] += trueRMS * noiseA[i];
    }
  }

  vector<float> a;
  for(int i=0; i<NSAMPLES; i++) a.push_back((float)amp[i]);
  return a;
} 

void test()
{
     
     
     double pedrms = 1;
     PulseChiSqSNNLS pulsefunc;
     PulseChiSqSNNLS pulsefuncSN;
     PulseChiSqSNNLS pulsefuncRef;
     
     // loop over ampSignal, true signal amplitude
     float ampSignalStart=0., ampSignalEnd=10., ampSignalStep=1.;
     for(float ampSignal=ampSignalStart; ampSignal<ampSignalEnd; ampSignal+=ampSignalStep){
	  digi signal(ampSignal);
	  float ampPUStart=0., ampPUEnd=10., ampPUStep=1.;
	  for(float ampPU=ampPUStart; ampPU<ampPUEnd; ampPU+=ampPUStep){
	       int bxPUStart=-10, bxPUEnd=10, bxPUStep=1;
	       for(int bxPU=bxPUStart; bxPU<bxPUEnd; bxPU+=bxPUStep){
		    
		    signal.addOOT(ampPU, bxPU);
		    vector<float> refAmp = signal.getDigis();
		    for(unsigned int i=0; i < NSAMPLES; i++){
			 amplitudes[i]=signal[i];
		    }
		    
		    
		    bool status = pulsefunc.DoFit(amplitudes,noisecor,pedrms,activeBX,fullpulse,fullpulsecov);
		    //  double chisq = pulsefunc.ChiSq();
		    PulseVector x = pulsefunc.X();
		    for( int i=0; i< NSAMPLES; ++i){
			 std::cout << i << "\t" << x[i] << "\t" << ampSignal << "\t" << ampPU << "\t" << bxPU << std::endl;;
		    }
	       }
	  }
     }
}

void saveHist()
{

  fout->cd();
//   h01->Write();
//   h02->Write();
//   h03->Write();
//   h04->Write();
//   h05->Write();
//   h06->Write();
//   h07->Write();
//   h11->Write();
//   h12->Write();

/*
  hMF_ampRef_vs_ampTrue->Write();
  hMF_Delta_ampRef_ampNoise->Write();
  hMF_relDelta_ampRef_amp_vs_BX->Write();
  hMF_relDelta_ampNoise_amp_vs_BX->Write();
  hMF_relDelta_ampTrue_amp_vs_BX->Write();
  hMF_Ratio_ampRef_amp_vs_BX->Write();
  hMF_Ratio_ampNoise_amp_vs_BX->Write();
  hMF_Ratio_ampTrue_amp_vs_BX->Write();

  hWg_ampRef_vs_ampTrue->Write();
  hWg_Delta_ampRef_ampNoise->Write();
  hWg_relDelta_ampRef_amp_vs_BX->Write();
  hWg_relDelta_ampNoise_amp_vs_BX->Write();
  hWg_relDelta_ampTrue_amp_vs_BX->Write();
  hWg_Ratio_ampRef_amp_vs_BX->Write();
  hWg_Ratio_ampNoise_amp_vs_BX->Write();
  hWg_Ratio_ampTrue_amp_vs_BX->Write();
*/


  hMF_ampRef_vs_BX->Write();
  hMF_ampNoise_vs_BX->Write();
  hMF_ampPU_vs_BX->Write();
  hWg_ampRef_vs_BX->Write();
  hWg_ampNoise_vs_BX->Write();
  hWg_ampPU_vs_BX->Write();



//   hMF_ampRef_ampTrue->Write();
//   hMF_ampRef_ampNoise->Write();

//   hMF_ampRef_amp->Write();
//   hMF_ampNoise_amp->Write();
//   hMF_ampTrue_amp->Write();
//   hMF_ampRef_amp_Ratio->Write();
//   hMF_ampNoise_amp_Ratio->Write();
//   hMF_ampTrue_amp_Ratio->Write();

//   hWg_ampRef_ampTrue->Write();
//   hWg_ampRef_ampNoise->Write();

//   hWg_ampRef_amp->Write();
//   hWg_ampNoise_amp->Write();
//   hWg_ampTrue_amp->Write();
//   hWg_ampRef_amp_Ratio->Write();
//   hWg_ampNoise_amp_Ratio->Write();
//   hWg_ampTrue_amp_Ratio->Write();
  

  fout->Close();
}



# ifndef __CINT__
int main()
{

  std::cout << " NSAMPLES = " << NSAMPLES << std::endl; 
  std::cout << " NToys = " << NTOYS << std::endl;
  std::cout << " nOOT = " << nOOT << std::endl;
  std::cout << " startOOT = " << startOOT << std::endl;
  std::cout << " ADC2GeVEB = " << ADC2GeVEB << std::endl;
  std::cout << " ADC2GeVEE = " << ADC2GeVEE << std::endl;

  init();
//  run(10);
  test();
  //    std::cout << " >>> main returned run " << std::endl;
  //    return 7;
  saveHist();
  return 0;
}
# endif
