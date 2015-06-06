#ifndef pulseshape_hh
#define pulseshape_hh
/**
   class pulseShape pulseShape include/pulseShape.hh
   
   pulse shape defined with alfa-beta function

*/
class pulseShape
{
public:
pulseShape(float tZero, float alpha=1.63669, float beta=1.44825):
     _tZero(tZero),
     _alpha(alpha),
     _beta(beta){
     
};
     
     double eval(double x, float tZero, float alpha, float beta);
     
     
     inline  float eval(double x){
	  return eval(x, _tZero, _alpha, _beta);
     }  
     
private:
     float  _tZero; ///< poistion of the maximum 
     float  _alpha, _beta;
};


#endif
