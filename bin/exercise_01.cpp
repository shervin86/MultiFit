#include <iostream>
#include "PulseChiSqSNNLS.h"
#include "MultiFitToyMC.h"

#include "Riostream.h"
#include "TF1.h"
#include "TProfile.h"
#include "TH2.h"
#include "TFile.h"
 
using namespace std;

double tZero = 0.0;
double alpha = 1.63669;
double beta  = 1.43412;

FullSampleVector fullpulse(FullSampleVector::Zero());
FullSampleMatrix fullpulsecov(FullSampleMatrix::Zero());
SampleMatrix noisecor(SampleMatrix::Zero());
BXVector activeBX;
SampleVector amplitudes(SampleVector::Zero());

TF1 *fps;

TFile *fout;
TH2D *h01;
TH2D *h02;
TH2D *h03;
TH2D *h04;
TH2D *h05;
TProfile *h10;
TProfile *h11;
TProfile *h12;
TProfile *h13;
TProfile *h14;
TProfile *h15;
TProfile *h16;
TProfile *h17;
TProfile *h18;
TProfile *h19;


void initHist()
{
  fout = new TFile("output.root","recreate");
  h01 = new TH2D("h01", "dA",     50, -1.0, 4.0, 1000, -20.0, 20.0);
  h02 = new TH2D("h02", "err",    50, -1.0, 4.0, 1000,  0.0, 20.0);
  h03 = new TH2D("h03", "dA/err", 50, -1.0, 4.0, 1000, -10.0, 10.0);
  h04 = new TH2D("h04", "chi2",   50, -1.0, 4.0, 1000,   0.0, 100.0);
  h05 = new TH2D("h05", "log(chi2)",   50, -1.0, 4.0, 1000,  -2., 3.);

  h10 = new TProfile("h10", "amp",   50, -1.0, 4.0, "S");
  h11 = new TProfile("h11", "amp",   50, -1.0, 4.0, "S");
  h12 = new TProfile("h12", "amp",   50, -1.0, 4.0, "S");
  h13 = new TProfile("h13", "amp",   50, -1.0, 4.0, "S");
  h14 = new TProfile("h14", "amp",   50, -1.0, 4.0, "S");
  h15 = new TProfile("h15", "amp",   50, -1.0, 4.0, "S");
  h16 = new TProfile("h16", "amp",   50, -1.0, 4.0, "S");
  h17 = new TProfile("h17", "amp",   50, -1.0, 4.0, "S");
  h18 = new TProfile("h18", "amp",   50, -1.0, 4.0, "S");
  h19 = new TProfile("h19", "amp",   50, -1.0, 4.0, "S");

}

void init()
{
  initHist();

  fps = new TF1("fps",funPulseShape, -10., 100., 3);
  fps->SetParameters(tZero+2.0, alpha, beta);
  double pulseShapeTemplate[12];
  for(int i=0; i<12; i++) pulseShapeTemplate[i] = fps->Eval((double)i);
  for(int i=0; i<12; i++) pulseShapeTemplate[i] /= pulseShapeTemplate[2];
  for (int i=0; i<12; ++i) fullpulse(i+7) = pulseShapeTemplate[i];

  for (int i=0; i<10; ++i) {
    for (int j=0; j<10; ++j) {
      if(i==j){
	noisecor(i,j) = 1;
      }else{
	noisecor(i,j) = 0;
      }
    }
  }

  int activeBXs[] = { -5, -4, -3, -2, -1,  0,  1,  2,  3,  4 };
  activeBX.resize(10);
  for (unsigned int ibx=0; ibx<10; ++ibx) {
    activeBX.coeffRef(ibx) = activeBXs[ibx];
  } 
}


vector<int> newHit(double aMax, double trueRMS, double pedestal)
{
  double amp[10];
  for(int i=0; i<10; i++) amp[i] = 0.;

  fps->SetParameters(tZero + 5.0, alpha, beta);
  for(int i=0; i<10; i++) amp[i] += aMax * fps->Eval((double)i);

  for(int i=0; i<10; i++){
    amp[i] += pedestal;
    amp[i] += rnd.Gaus(0., trueRMS);
  }

  vector<int> a;
  for(int i=0; i<10; i++) a.push_back((int)amp[i]);
  return a;
}
 

void run()
{
  int ntot = 100000;
  for(int ievt=0; ievt<ntot; ++ievt){

    double ampTrue = pow(10, -1.0 + 5.0 * rnd.Rndm());
    double pedTrue = 180.0 +  40.0 * rnd.Rndm();
    double rmsTrue = 2.0;
    vector<int> amp = newHit(ampTrue, rmsTrue, pedTrue);

    for(int i=0; i<10; i++){
      amplitudes[i] = amp[i] - (pedTrue - 0.5);
    }

    double pedval = 0.;
    double pedrms = sqrt(rmsTrue * rmsTrue + 1./12.) ;  
    PulseChiSqSNNLS pulsefunc;

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

    double logA = log10(ampTrue);
    h01->Fill( logA, aMax - ampTrue);
    h02->Fill( logA, aErr );
    if(aErr>0) h03->Fill( logA, (aMax - ampTrue)/aErr );
    h04->Fill( logA, chisq );
    if(chisq>0) h05->Fill( logA, log10(chisq) );
    
    h10->Fill( logA, amplitudes[0] );
    h11->Fill( logA, amplitudes[1] );
    h12->Fill( logA, amplitudes[2] );
    h13->Fill( logA, amplitudes[3] );
    h14->Fill( logA, amplitudes[4] );
    h15->Fill( logA, amplitudes[5] );
    h16->Fill( logA, amplitudes[6] );
    h17->Fill( logA, amplitudes[7] );
    h18->Fill( logA, amplitudes[8] );
    h19->Fill( logA, amplitudes[9] );

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
  h10->Write();
  h11->Write();
  h12->Write();
  h13->Write();
  h14->Write();
  h15->Write();
  h16->Write();
  h17->Write();
  h18->Write();
  h19->Write();
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
