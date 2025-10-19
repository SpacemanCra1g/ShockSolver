#include "../include/DomainClass.hpp"

int Domain::FixedVx(double *Uin, double *Uout, const int i, const double vx){

  double D = Uin[Tidx(DENS, i)];
  double Sx = Uin[Tidx(MOMX, i)];
  double Sy = Uin[Tidx(MOMY, i)];
  double E = Uin[Tidx(ENER, i)];

  double schr = Sx/vx;
  double vy = Sy/schr;
  double p = schr - E;
  double rho = D*std::sqrt(1.0 - vx*vx - (Sy/schr)*(Sy/schr));

  Uout[Tidx(DENS,i)] = rho;
  Uout[Tidx(VELX,i)] = vx;
  Uout[Tidx(VELY,i)] = vy;
  Uout[Tidx(PRES,i)] = p;


  return 0;
}
