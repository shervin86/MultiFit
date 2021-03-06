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

#define NSAMPLES 10
#define NToys 1000
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



void run()
{
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
  PulseChiSqSNNLS pulsefuncSN;
  PulseChiSqSNNLS pulsefuncRef;

  //FIXME need to decide nToys
  for(int iToy=0; iToy<NToys; ++iToy){
    if(iToy % 200 == 0) std::cout << " >>> iToy = " << iToy << std::endl;

    //NB start witn fix amplitude, no ped, no noise

    //    double ampSignal = pow( 10., -1.0 + 5.0 * rnd.Rndm() );
    double ampSignal = 1.;
    double ampRef = ampSignal;

    //FIX why true and not real?
    vector<float> nullPU = newDigiPU(0., 0.);
    //    vector<int> refAmp = newDigi(0., 0., pedTrue, ampSignal, 0);
    vector<float> refAmp = newDigi(0., 0., 0., ampSignal, nullPU, 0);
    //    vector<int> ampNoise = newDigi(0., rmsTrue, pedTrue, ampSignal, 0);
    vector<float> ampNoise = newDigi(0., 0., 0., ampSignal, nullPU, 0);

    //FIXME cambia loop su OOT PU -> estendendolo da +2 <> -13 e poi +3 <> -13   e +4 <> -13
    for(int iOOT=startOOT; iOOT<= (stopOOT); ++iOOT){
      //      if(iToy % 200 == 0) std::cout << " >>> iOOT = " << iOOT << std::endl;

      vector<float> amp;
      vector<float> ootPU;
      amp.clear();
      ootPU.clear();

      //adding PU in window
      for(int iOOTinit=iOOT; iOOTinit<= (stopOOT); ++iOOTinit){
	//	double ampPUIT   = pow( 10., -1.0 + 5.0 * rnd.Rndm() );
	double ampPUIT = 1.;
	double ampTrue = ampSignal + ampPUIT;  // useless

	vector<float> ootPU_loc = newDigiPU(ampPUIT, iOOT);
	for(int ii=0; ii<NSAMPLES; ++ii){
	  if(iOOTinit == iOOT) ootPU.push_back(ootPU_loc.at(ii));
	  else ootPU.at(ii) += ootPU_loc.at(ii);
	}
	//	amp = newDigi(ampPUIT, rmsTrue, pedTrue, ampSignal, iOOTinit);
	//amp = newDigi(0., 0., 0., ampSignal, ootPU, iOOTinit);
	amp = newDigi(0., 0., 0., 0., ootPU, iOOTinit);
      }

    //FIXME why true and not real?
    for(int i=0; i<NSAMPLES; i++){
      //      amplitudes[i] = amp[i] - (pedTrue - 0.5);
      amplitudes[i] = amp[i];
      //      amplitudesNoise[i] = ampNoise[i] - (pedTrue - 0.5);
      amplitudesNoise[i] = ampNoise[i];
      //      amplitudesRef[i] = refAmp[i] - (pedTrue - 0.5);
      amplitudesRef[i] = refAmp[i];
    }
    
    if(iOOT == -2 && iToy == 1){
    TH1F* h_amp = new TH1F("h_amp", "", 30, 0., 30.);
    TH1F* h_ampN = new TH1F("h_ampN", "", 30, 0., 30.);
    TH1F* h_ampR = new TH1F("h_ampR", "", 30, 0., 30.);
    TH1F* h_fullP = new TH1F("h_fullP", "", 30, 0., 30.);
    
    for(int i=0; i<NSAMPLES; i++){
      std::cout << " >>> amplitudes[i] = " << amplitudes[i] << std::endl;
      h_amp->Fill(i, amplitudes[i]);
      std::cout << " >>> amplitudesNoise[i] = " << amplitudesNoise[i] << std::endl;
      h_ampN->Fill(i, amplitudesNoise[i]);
      std::cout << " >>> amplitudesRef[i] = " << amplitudesRef[i] << std::endl;
      h_ampR->Fill(i, amplitudesRef[i]);
    }
    for(int i=0; i<19; i++){
      std::cout << " >>> fullpulse[i] = " << fullpulse[i] << std::endl;
      h_fullP->Fill(i, fullpulse[i]);
    }

    std::cout << " >>> while looping iOOT = " << iOOT << std::endl;

    TFile out("histos.root", "recreate");
    out.cd();
    h_amp->Write("h_amp");
    h_ampN->Write("h_ampN");
    h_ampR->Write("h_ampR");
    h_fullP->Write("h_fullP");
    out.Close();
    std::cout << " >>> before return " << std::endl;
    //    return;
    }

    //    if(iToy % 200 == 0) std::cout << " >>> start fitting " << std::endl;
    // fit signal + noise + PU
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

    //    if(iToy % 200 == 0) std::cout << " >>> start fitting " << std::endl;
    // fit signal + noise 
    //  pulsefunc.disableErrorCalculation();
    bool statusSN = pulsefuncSN.DoFit(amplitudesNoise, noisecor,pedrms,activeBX,fullpulse,fullpulsecov);
    double chisqSN = pulsefuncSN.ChiSq();
  
    ipulseintime = 0;
    for (unsigned int ipulse=0; ipulse<pulsefuncSN.BXs().rows(); ++ipulse) {
      if (pulsefuncSN.BXs().coeff(ipulse)==0) {
	ipulseintime = ipulse;
	break;
      }
    }
    double aMaxSN = statusSN ? pulsefuncSN.X()[ipulseintime] : 0.;
    double aErrSN = statusSN ? pulsefuncSN.Errors()[ipulseintime] : 0.;

    double aWgtSN = 0.;
    for(int i=0; i<NSAMPLES; i++){
      aWgtSN += weights[i] * amplitudesNoise[i];
    }

    //    if(iToy % 200 == 0) std::cout << " >>> start fitting " << std::endl;
    // fit signal
    //  pulsefunc.disableErrorCalculation();
    bool statusR = pulsefuncRef.DoFit(amplitudesRef, noisecor,pedrms,activeBX,fullpulse,fullpulsecov);
    double chisqR = pulsefuncRef.ChiSq();
  
    ipulseintime = 0;
    for (unsigned int ipulse=0; ipulse<pulsefuncRef.BXs().rows(); ++ipulse) {
      if (pulsefuncRef.BXs().coeff(ipulse)==0) {
	ipulseintime = ipulse;
	break;
      }
    }
    double aMaxR = statusR ? pulsefuncRef.X()[ipulseintime] : 0.;
    double aErrR = statusR ? pulsefuncRef.Errors()[ipulseintime] : 0.;

    double aWgtR = 0.;
    for(int i=0; i<NSAMPLES; i++){
      aWgtR += weights[i] * amplitudesRef[i];
    }



    //    if(iToy % 200 == 0) std::cout << " >>> start filling " << std::endl;

    /*
    double logA = log10(ampRef * 0.069);

    hMF_ampRef_vs_ampTrue->Fill(ampRef, aMaxR);
    hMF_Delta_ampRef_ampNoise->Fill((aMaxR - aMaxSN) * 0.069);

    //    if(aMaxSN+aMax == 0) std::cout << " >>>> ahhhhhh " << std::endl;

    if(aMaxR != 0.)    hMF_relDelta_ampRef_amp_vs_BX->Fill(iOOT, (aMax - aMaxR) / (aMaxR));
    if(aMaxSN != 0.)    hMF_relDelta_ampNoise_amp_vs_BX->Fill(iOOT, (aMax - aMaxSN) / (aMaxSN));
    if(ampRef != 0.)    hMF_relDelta_ampTrue_amp_vs_BX->Fill(iOOT, (aMax - ampRef) / (ampRef));

    if(aMaxR != 0.) hMF_Ratio_ampRef_amp_vs_BX->Fill(iOOT, (aMax/aMaxR));

    //analogamente... sara' questo?  plot e' l'unico ok
    if(aMaxSN != 0.) hMF_Ratio_ampNoise_amp_vs_BX->Fill(iOOT, (aMax/aMaxSN));

    if(ampRef != 0.) hMF_Ratio_ampTrue_amp_vs_BX->Fill(iOOT, (aMax/ampRef));

    //    hWg_ampRef_ampTrue->Fill(logA, aWgtR * 0.069);
    hWg_ampRef_vs_ampTrue->Fill(ampRef, aWgtR);
    hWg_Delta_ampRef_ampNoise->Fill((aWgtR - aWgtSN) * 0.069);

    if(aWgtR != 0.)    hWg_relDelta_ampRef_amp_vs_BX->Fill(iOOT, (aWgt - aWgtR) / (aWgtR));
    if(aWgtSN != 0.)    hWg_relDelta_ampNoise_amp_vs_BX->Fill(iOOT, (aWgt - aWgtSN) / (aWgtSN));
    if(ampRef != 0.)    hWg_relDelta_ampTrue_amp_vs_BX->Fill(iOOT, (aWgt - ampRef) / (ampRef));

    if(aWgtR != 0) hWg_Ratio_ampRef_amp_vs_BX->Fill(iOOT, aWgtR);

    //questo plot e' l'unico ok
    if(aWgtSN != 0) hWg_Ratio_ampNoise_amp_vs_BX->Fill(iOOT, (aWgt/aWgtSN));
    //    std::cout << " >>> iOOT = " << iOOT << " aWgt/aWgtSN = " << aWgt/aWgtSN << std::endl;

    //fixme
    if(ampRef != 0) hWg_Ratio_ampTrue_amp_vs_BX->Fill(iOOT, aWgt);
    */


    hMF_ampRef_vs_BX->Fill(iOOT, aMaxR);
    hMF_ampNoise_vs_BX->Fill(iOOT, aMaxSN);
    hMF_ampPU_vs_BX->Fill(iOOT, aMax);


    hWg_ampRef_vs_BX->Fill(iOOT, aWgtR);
    hWg_ampNoise_vs_BX->Fill(iOOT, aWgtSN);
    hWg_ampPU_vs_BX->Fill(iOOT, aWgt);


    
    /*  
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
    */

    } // looping over PU events
  }//loop over toys

  std::cout << " >>> end call " << std::endl;
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
  std::cout << " NToys = " << NToys << std::endl;
  std::cout << " nOOT = " << nOOT << std::endl;
  std::cout << " startOOT = " << startOOT << std::endl;
  std::cout << " ADC2GeVEB = " << ADC2GeVEB << std::endl;
  std::cout << " ADC2GeVEE = " << ADC2GeVEE << std::endl;

  init();
  run();
  //    std::cout << " >>> main returned run " << std::endl;
  //    return 7;
  saveHist();
  return 0;
}
# endif
