/*
 Harm Schoorlemmer May 13, 2014.  harmscho@phys.hawaii.edu
 
 -=Library to calculate atmospheric properties for
 extensive air shower.=-
 
 Currently only a Linsley density
 profile is implemented, based on the parameterization 
 given in the Aires Manual (http://www2.fisica.unlp.edu.ar/auger/aires/eg_Aires.html). 
 By default distance and altitude are calculated for a 
 curved Earth geometry, special functions are provided to
 compare those to a flat surface geometry. 
 For the refractivity model there are several options, the default corresponds to the model which is implemented in 
 ZHAireS.
 
 Do not distribute this library without my premission and use it at own risk.
 
*/
#include <math.h>
//parsec in meters
//all in SI-base units

const double cRad2Deg = 180./ M_PI;
const double cR_earth = 6371e3;
const double cc = 3e8;


//conversions to seconds.
const double cMin = 60;
const double cHour = 60 * cMin;
const double cDay = 24 * cHour;
const double cYear = 365.24 * cDay;

const double cKilo = 1e3;
const double cMega = 1e6;
const double cMili = 1e-3;
const double cMicro = 1e-6;

void SayHi();

enum DensityModel {
  eLinsley
};

class Atmosphere {
public:

  
  Atmosphere(DensityModel m);

  void CurrentModel();
  
  double GetDensityAtHeight(double h);
  double GetVerticalDepthAtHeight(double h);
  double VerticalHeight(double rho, double z, double z_gr);
  double GetSlantDepth(double zenith,double r,double z_gr,bool useCurved = true);
  double GetDistanceToSlantDepthFlatEarth(double zenith,double X, double z_gr);
  double GetDistanceToSlantDepth(double zenith,double X, double z_gr);
  double RefractivityAtHeight(double h); // above sea level
  double RefractivityAtHeightModified(double h); // above sea level
  double RefractivityAtHeightSouthPole(double h); // above sea level
  double GetSlantDepthEmission(double rmax,double zenith,double z_gr);//rmax is the lateral distance where the maximum emission occurs
  
private:
  double GetVerticalDepthAtHeightLinsley(double h /* m */);//returns g / cm^2
  double GetDensityAtHeightLinsley(double h /* m */);//returns g / cm^3
  DensityModel fModel;

};






