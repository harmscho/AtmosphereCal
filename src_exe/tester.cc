#include "Atmosphere.h"
#include <iostream>
#include <fstream>
using namespace std;

int main() {

  Atmosphere atm(eLinsley);
  
  
  double d2r = M_PI / 180;//degree to rad
  double zGround = 1400; // 50 m
  double zenithAngle = 0 * d2r;//zenith angle
  double slantDepth = 448; // g / cm^2
  
  cout << "\n======Example=====" << endl;
  cout << "Using a Linsley atmosphere,\naltitude of ground = 50 m a. s. l,\nZenith angle of 70 deg.,\nSlant depth of Xmax is 650 g/cm^2 " << endl;
  cout << "\n------Results-----" << endl;

  double dXmax = atm.GetDistanceToSlantDepth(zenithAngle,slantDepth,zGround);
  cout << "Distance from impact point on ground to Xmax is: " << dXmax/1e3  << " [km]" << endl;

  double altitudeXmax = atm.VerticalHeight(dXmax * sin(zenithAngle), dXmax * cos(zenithAngle),zGround);
  cout  << "Xmax is at an altitude of " <<altitudeXmax/1e3 << "[km] above sea level." << endl;
  
  double eta = atm.RefractivityAtHeight(altitudeXmax);
  double psi_C = acos(1./(1+ eta));
  
  cout << "The refractivy at Xmax is " << eta << "," <<endl;
  cout << "resulting in a Cherenkov angle of " << psi_C/d2r << " [deg]\n\n" <<  endl;
  
  return 0;
}