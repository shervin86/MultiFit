#ifndef digi_hh
#define digi_hh

#include <vector>
#include <ostream> 

#ifndef MAX_SAMPLES
#define MAX_SAMPLES 10
#endif

class digi {
public:
     /// default constructor
     digi():
	  _signalTime(5.0),
	  _pedestal(0.){}; 

     digi(double aMax):
	  _signalTime(5.0),
	  _pedestal(0.){
	  setSignal(aMax);
     }; 


     void setSignal(double aMax); ///< set the digis to the values given by pulseShape with max = aMax and at the 6th sample
     void setOOT(double aMax, double iBX);
     void addSignal(double aMax); ///< add the values given by pulseShape with max = aMax and at the 6th sample to the digis 
     void addOOT(double aMax, double iBX);
     inline void setPedestal(double pedestal){
	  _pedestal=pedestal;
     }

     /// return a vector with the values saved in digis
     inline std::vector<float> getDigis(void){
	  std::vector<float> a;
	  for(int i=0; i<MAX_SAMPLES; i++) a.push_back((float) digis[i]+noise[i]+_pedestal);
	  return a;
     }
     
     inline float operator[](unsigned int i) const{
	  return digis[i];
     }
     inline friend digi operator +(digi& lhs, digi& rhs){
	  digi newdigi;
	  for(unsigned int i=0; i< MAX_SAMPLES; ++i){
	       newdigi.digis[i]=rhs.digis[i]+lhs.digis[i];
	  }
	  return newdigi;
     }

     inline friend void operator +=(digi& lhs, digi &rhs){
	  for(unsigned int i=0; i< MAX_SAMPLES; ++i){
	       lhs.digis[i]+=rhs.digis[i];
	  }
     }


     inline friend std::ostream& operator <<(std::ostream& os, digi& m){
	  for(unsigned int i=0; i< MAX_SAMPLES; ++i){
	       os << "digi:\t" << m.digis[i] << "\n";
	  }
	  return os;
     }	  

private:
     float _signalTime;
     float digis[MAX_SAMPLES];
     float noise[MAX_SAMPLES];
     float _pedestal;
     bool hasPedestal;
     bool hasNoise;
     


};

#endif
