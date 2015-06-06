#include "digi.hh"
#include "pulseShape.hh"

void digi::setSignal(double aMax){
  //  building pulse shape to be fitted
  
  //add signal
  pulseShape ps(_signalTime);
 
  for(unsigned int i=0; i<MAX_SAMPLES; i++)     digis[i] = aMax * ps.eval((double)i);
    return;
}

void digi::setOOT(double aMax, double tZero){
     //  building pulse shape to be fitted
     
     //add signal
     pulseShape ps(tZero);
     
     for(unsigned int i=0; i<MAX_SAMPLES; i++)     digis[i] = aMax * ps.eval((double)i);
     return;
}


void digi::addOOT(double aMax, double tZero){
     pulseShape ps(tZero);
     for(unsigned int i=0; i<MAX_SAMPLES; i++)     digis[i] += aMax * ps.eval((double)i);
     return;
}

void digi::addSignal(double aMax){
     pulseShape ps(_signalTime);
     for(unsigned int i=0; i<MAX_SAMPLES; i++)     digis[i] += aMax * ps.eval((double)i);
     return;
}


