#include "Atmosphere.h"
#include <iostream>
#include <fstream>
using namespace std;

void SayHi() {
  cout << "Hello!" << endl;

}

Atmosphere::Atmosphere(DensityModel m) {
   fModel = m;
}

void Atmosphere::CurrentModel()
{
  if (fModel == eLinsley) {
    cout << "Atmosphere is a Linsley model" << endl;
  } else cout << "Model unknown" << endl;
}

double Atmosphere::GetDensityAtHeight(double h)
{
  if (fModel == eLinsley) {
    return GetDensityAtHeightLinsley(h);
  } else cout << "Model unknown" << endl;
  return 0;
}

double Atmosphere::GetVerticalDepthAtHeight(double h)
{
  if (fModel == eLinsley) {
    return GetVerticalDepthAtHeightLinsley(h);
  } else cout << "Model unknown" << endl;
  
  return 0;
}

double Atmosphere::GetDensityAtHeightLinsley(double h)
{
  return -(GetVerticalDepthAtHeightLinsley(h + 0.5) - GetVerticalDepthAtHeightLinsley(h - 0.5))/100;
}


double Atmosphere::GetVerticalDepthAtHeightLinsley(double h)
{
  
  const double a[6] = {-186.5562,-94.9199,0.61289,0,1.0469564e-10,2.614186472e-10};//g / cm^2
  const double b[6] = {1222.6562,1144.9069,1305.5948,540.1778,540.17833,1.};////g / cm^2
  const double c[6] = {9941.8638,8781.5355,6361.4304,7721.7016,7721.70097,1.64444076e15}; // m
  const double lim [6] = {4.0e3,10.0e3,40.0e3,100.0e3,250e3,429.89e3}; // m
	
  if (h < lim[0]) {
		//cout << 0 << endl;
		return a[0] + (b[0] * exp(-h/c[0]));
	} else if(h >= lim[0] && h < lim[1]) {
		//cout << 1 << endl;
		return a[1] + (b[1] * exp(-h/c[1]));
	} else if(h >= lim[1] && h < lim[2]) {
		//cout << 2 << endl;
		return a[2] + (b[2] * exp(-h/c[2]));
	} else if(h >= lim[2] && h < lim[3]) {
		return a[3] + (b[3] * exp(-h/c[3]));
	} else if(h >= lim[3] && h < lim[4]) {
		return a[4] + (b[4] * exp(-h/c[4]));
	} else if(h >= lim[4] && h < lim[5]) {
		return a[5] - (b[5] * h /c[5]);
		cout << h << endl;
	} else if(h >= lim[5]) {
		return 0;
	} else {
		cout << "ERROR WTF is a h of " << h << endl;
		return 0;
	}
  return 0;
}



double Atmosphere::VerticalHeight(double rho, double z, double z_gr)
{
  double z_c = z + z_gr;
  double C = -(pow(cR_earth + z_c,2) + pow(rho,2) - pow(cR_earth,2));
  double B = 2 * cR_earth;
  
  double z_v_min = (-B - sqrt(pow(B,2) - (4 * C))) / 2;
  double z_v_plus = (-B + sqrt(pow(B,2) - (4 * C))) / 2;
  if (z_v_min > 0) cout << "z_v_min " << z_v_min << "\t z_v_plus " << z_v_plus << endl;
  return z_v_plus;
}

double Atmosphere::GetSlantDepth(double zenith,double r,double z_gr,bool useCurved)
{
 
  double dr = 5.;
  double N = 0.2e6;
  double sum = 0;
  double rho,z,z_v;
  double _r;
  for (int i = 0; i < N; i++) {
    _r = r + i * dr;
    rho = sin(zenith) * _r;
    z = cos(zenith) * _r;
    if (useCurved) {
      z_v = VerticalHeight(rho,z,z_gr);
    } else {
      z_v = z + z_gr;
    }
    
    sum += GetDensityAtHeight(z_v) * dr;
  }
  return sum * 100;//g / cm^2
}

double Atmosphere::GetDistanceToSlantDepth(double zenith,double X, double z_gr) {

  double _X = 1e30;
  double delta = 1e6;
  double _r = 0;
  int cnt = 0;
  
  if (X > GetSlantDepth(zenith, 0, z_gr)) {
    cout << "WARNING: requested slant depth ( "<< X <<" g / cm^2) is larger than slant depth at ground level ( " <<GetSlantDepth(zenith,0, z_gr)<< " g / cm^2)"  << endl;
    cout << "Setting distance to Slant Depth to -1" << endl;
    return -1;
  }
  
  
  
  while (fabs(_X - X) > 1) {
    delta /= 2;
    if (_X > X) {
      _r = _r + delta;
    } else {
      _r = _r - delta;
    }
    _X = GetSlantDepth(zenith, _r, z_gr);
    cnt++;
    if (cnt > 20) {
      cout << "WARNING: GetDistanceToSlantDepth did not converged!!" << endl;
      return -1;
    }
  }
  
  return _r;
}

double Atmosphere::GetDistanceToSlantDepthFlatEarth(double zenith,double X, double z_gr) {
  
  double _X = 1e30;
  double delta = 1e6;
  double _r = 0;
  int cnt = 0;
  while (fabs(_X - X) > 1) {
    delta /= 2;
    if (_X > X) {
      _r = _r + delta;
    } else {
      _r = _r - delta;
    }
    _X = GetVerticalDepthAtHeight(_r * cos(zenith) + z_gr)/cos(zenith);
    cnt++;
    //cout << _X << "\t" << X << endl;
  }
  
  return _r;
}

double Atmosphere::RefractivityAtHeight(double h)
{
  double kr = -0.1218;//km^1
  double n0 = 325e-6;
  h /= 1e3;
  return n0 * exp(kr * h);


}

double Atmosphere::RefractivityAtHeightModified(double h)
{
  double kr = -0.094;//km^1
  double n0 = 340e-6;
  h /= 1e3;
  return n0 * exp(kr * h);  
}


double Atmosphere::RefractivityAtHeightSouthPole(double h)
{
  double kr = -0.1369;//km^1
  double n0 = 327.4e-6;
  h /= 1e3;
  return n0 * exp(kr * h);
}


double Atmosphere::GetSlantDepthEmission(double rmax,double zenith,double z_gr)
{
  double X,z_v;
  double rho,z;
  double _rmax = 1e80;
  double rLong;
  double nr,psi;
  ofstream ofs("rl_rt.txt");
  //There is a maximum radius at which the radiation can end up
  double max = 0;
  double rLongMax;
  for (int i = 0; i < 100000; i++) {
    rLong = i * 50;
    rho = sin(zenith) * rLong;
    z = cos(zenith) * rLong;
    z_v = VerticalHeight(rho, z, z_gr);
    nr = RefractivityAtHeightModified(z_v);
    psi = acos(1./(1. + nr));
    _rmax = tan(psi) * rLong;
     ofs << rLong << "\t" << _rmax << endl;
    if (_rmax > max) {
     max = _rmax;
      rLongMax = rLong;
    }
  }

  if ( rmax > max ) {
    cout << "The refractive angle cannot explain a maximum at " << rmax << endl;
    cout << " Returing -1" << endl;
    return -1;
  }

  cout << "Maximum cherenkov radius: " << max << " m, at logitudinal distance: " << rLongMax << endl;
  //finding the emmison point
  int cnt = 0;
  _rmax = max;
  double delta = rLongMax;
  rLong = delta;

  while (fabs(_rmax - rmax) > 1 && cnt < 20) {
    delta /= 2;
    if (_rmax < rmax) {
      rLong = rLong + delta;
    } else {
      rLong = rLong - delta;
    }
    
    rho = sin(zenith) * rLong;
    z = cos(zenith) * rLong;
    z_v = VerticalHeight(rho, z, z_gr);
    nr = RefractivityAtHeight(z_v);
    psi = acos(1./(1. + nr));
    _rmax = tan(psi) * rLong;
    cnt++;

    cout << cnt << "\t" << _rmax << "\t" << z_v/1e3 << "\t" << rLong/1e3 << endl;
  }

  if (cnt >= 20) {
    cout << "WARNING to many iterations needed" <<  endl;
  }
  
  
  X = GetSlantDepth(zenith,rLong,z_gr);
  return X;

  

}


