#include <iostream>
#include "PulseChiSqSNNLS.h"
#include "MultiFitToyMC.h"

#include "Riostream.h"
#include "TTree.h"
#include "TF1.h"
#include "TProfile.h"
#include "TH2.h"
#include "TFile.h"
 
using namespace std;

#define NSAMPLES 10

double tZero = 0.0;
double alpha = 1.63669;
double beta  = 1.43412;
double weights[NSAMPLES] = {
  -0.409783, -0.409783, -0.409783, 0.0, 0.258171,
  0.318941, 0.283887, 0.216833, 0.151517, 0.0 
};

FullSampleVector fullpulse(FullSampleVector::Zero());
FullSampleMatrix fullpulsecov(FullSampleMatrix::Zero());
SampleMatrix noisecor(SampleMatrix::Zero());
BXVector activeBX; // definito dove?
SampleVector amplitudes(SampleVector::Zero());

double EECorrNoiseMatrixG12[] = {
  1.00000, 0.71373, 0.44825, 0.30152, 0.21609,
  0.14786, 0.11772, 0.10165, 0.09465, 0.08098 
};
double LGen[NSAMPLES][NSAMPLES];

TF1 *fps;

TFile *fout;
TH2D *h01;
TH2D *h02;
TH2D *h03;
TH2D *h04;
TH2D *h05;
TH2D *h06;
TH2D *h07;

TH1D *h11;
TH1D *h12;

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
  fout = new TFile("output.root","recreate");
  h01 = new TH2D("h01", "dA",     50, -2.0, 3.0, 1000, -2.0, 2.0);
  h02 = new TH2D("h02", "err",    50, -2.0, 3.0, 1000,  0.0, 2.0);
  h03 = new TH2D("h03", "dA/err", 50, -2.0, 3.0, 1000, -10.0, 10.0);
  h04 = new TH2D("h04", "chi2",   50, -2.0, 3.0, 1000,   0.0, 100.0);
  h05 = new TH2D("h05", "log(chi2)",   50, -2.0, 3.0, 1000,  -2., 6.);
  h06 = new TH2D("h06", "dA",     50, -2.0, 3.0, 1000, -10.0, 10.0);
  h07 = new TH2D("h07", "dA",     50, -2.0, 3.0, 1000, -10.0, 10.0);

  h11 = new TH1D("h11", "IT PU",   1000,  0.0, 50.0);
  h12 = new TH1D("h12", "IT - OOT PU", 1000,  -20.0, 20.0);

}


void init()
{
  initHist();
  NoiseCorrInit();

  fps = new TF1("fps",funPulseShape, -10., 100., 3);
  fps->SetParameters(tZero+2.0, alpha, beta);
  double pulseShapeTemplate[12];
//  fps->Draw();
  //fps->SaveAs("pulseShape.root");

  for(int i=0; i<12; i++) pulseShapeTemplate[i] = fps->Eval((double)i);
  for(int i=0; i<12; i++) pulseShapeTemplate[i] /= pulseShapeTemplate[2]; //normalized to maximum
  for (int i=0; i<12; ++i) fullpulse(i+7) = pulseShapeTemplate[i];

  unsigned int num_activeBXs=NSAMPLES;

  for (int i=0; i<NSAMPLES; ++i) {
    for (int j=0; j<NSAMPLES; ++j) {
      int vidx = std::abs(j-i);
      noisecor(i,j) = EECorrNoiseMatrixG12[vidx];
    }
  }

  int activeBXs[] = { -5, -4, -3, -2, -1,  0,  1,  2,  3,  4 };
  activeBX.resize(num_activeBXs); 

  for (unsigned int ibx=0; ibx<num_activeBXs; ++ibx) {
    activeBX.coeffRef(ibx) = activeBXs[ibx];
  } 
  //  activeBX.resize(1);
  //  activeBX.coeffRef(0) = 0;
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
 


vector<int> newDigi( const double *e, double trueRMS, double pedestal, double aMax, int tMax=0 ) // tMax=0 in-time signal
{
  double amp[NSAMPLES];
  //for(int i=0; i<NSAMPLES; i++) amp[i] = e[i] / 0.069;

  fps->SetParameters(tMax + 5.0, alpha, beta);
  for(int i=0; i<NSAMPLES; i++) amp[i] += aMax * fps->Eval((double)i);


  if(pedestal>0.){
  for(int i=0; i<NSAMPLES; i++){
    amp[i] += pedestal;
  }
  }

  if(trueRMS>0.){
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



void run()
{
  TFile *fin1 = new TFile("minbiasEE_digi_eta28_pu50.root");
  TTree *tr = (TTree*)fin1->Get("PU");
  double emb[25];
  double ene[NSAMPLES];
  tr->SetBranchAddress("emb",emb);
  tr->SetBranchAddress("ene",ene);


  double pedTrue = 200.0;
  double rmsTrue = 2.0;

  for (int i=0; i<NSAMPLES; ++i) {
    for (int j=0; j<NSAMPLES; ++j) {
      int vidx = std::abs(j-i);
      if(vidx==0){
	noisecor(i,j) = 1.0;
      }else{
	noisecor(i,j) = EECorrNoiseMatrixG12[vidx] * rmsTrue * rmsTrue / ( rmsTrue * rmsTrue + 1./12. );
      }
    }
  }

  double pedval = 0.;
  double pedrms = sqrt(rmsTrue * rmsTrue + 1./12.) ;  
  //    double pedrms = 1.0;
  PulseChiSqSNNLS pulsefunc;

  int ntot = tr->GetEntries();
  for(int ievt=0; ievt<ntot; ++ievt){

    tr->GetEntry(ievt);

    h11->Fill( emb[20] );
    double sum0 = 0;
    double sum1 = 0;
    for(int i=0; i<25; i++){
      sum0 += 1;
      sum1 += emb[i];
    }
    sum0 -= 1.;
    sum1 -= emb[20];
    h12->Fill( emb[20] - sum1 / sum0 );

    double ampPUIT   = emb[20] / 0.069;
    double ampSignal = pow( 10., -1.0 + 5.0 * rnd.Rndm() );
    double ampTrue = ampSignal + ampPUIT;
    vector<int> amp = newDigi(ene, rmsTrue, pedTrue, ampSignal);

    
    for(unsigned int iOOT=0; iOOT< nOOT; ++iOOT){
	 
	 vector<int> oot_amp = newDigi(ene, rmsTrue, pedTrue, ampSignal);
    for(int i=0; i<NSAMPLES; i++){
      amplitudes[i] = amp[i] - (pedTrue - 0.5);
    }

    //  pulsefunc.disableErrorCalculation();
    bool status = pulsefunc.DoFit(amplitudes,noisecor,pedrms,activeBX,fullpulse,fullpulsecov);
    double chisq = pulsefunc.ChiSq();
  
    unsigned int ipulseintime = 0;
    for (unsigned int ipulse=0; ipulse<pulsefunc.BXs().rows(); ++ipulse) {
      if (pulsefunc.BXs().coeff(ipulse)==0) {
	ipulseintime = ipulse;
	break;
      }
    }
    double aMax = status ? pulsefunc.X()[ipulseintime] : 0.;
    double aErr = status ? pulsefunc.Errors()[ipulseintime] : 0.;

    double aWgt = 0.;
    for(int i=0; i<NSAMPLES; i++){
      aWgt += weights[i] * amplitudes[i];
    }
    
    double logA = log10(ampTrue * 0.069);
    h01->Fill( logA, (aMax - ampTrue) * 0.069 );
    h02->Fill( logA, aErr * 0.069 );
    if(aErr>0){
      h03->Fill( logA, (aMax - ampTrue)/aErr );
    }else{
      h03->Fill( logA, -19.0 );
    }
    h04->Fill( logA, chisq );
    if(chisq>0){
      double tmp = log10(chisq);
      if(tmp<-2) tmp = -1.9;
      if(tmp> 6) tmp = 5.9;
      h05->Fill( logA, tmp);
    }
    h06->Fill( logA, (aWgt - ampTrue) * 0.069 );
    h07->Fill( logA, (aWgt - ampSignal) * 0.069 );
  }  
}

void saveHist()
{

  fout->cd();
  h01->Write();
  h02->Write();
  h03->Write();
  h04->Write();
  h05->Write();
  h06->Write();
  h07->Write();
  h11->Write();
  h12->Write();
  fout->Close();
}



# ifndef __CINT__
int main()
{
  init();
  run();
  saveHist();
  return 0;
}
# endif
