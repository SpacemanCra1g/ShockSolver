#ifndef FLUXCLASS_H_
#define FLUXCLASS_H_

#include "../include/Parameters.h"
#include <cmath>
//  This is maybe a mistake, but I'm going to define a Flux container to
//  hold all of flux data objects and methods

class FluxClass {

public:
  double *Cons;
  double *dt;

  double FluxDir[NDIMS * 2][nqp][NumVar][xDim * yDim];

  double Flux[nqp][NumVar][xDim * yDim];

  void Fluxinit(double *COns, double *DT) {
    Cons = COns;
    dt = DT;
  }

  void FOG() {

    for (int qp = 0; qp < nqp; ++qp) {
      for (int var = 0; var < NumVar; ++var) {
        for (int xdir = 0; xdir < xDim; ++xdir) {
          for (int ydir = 0; ydir < yDim; ++ydir) {
            if (xdir > 0) {
              FluxDir[Left][qp][var][idx(xdir, ydir)] =
                  Cons[Tidx(var, xdir, ydir)];
            }

            if (xdir < xDim - 1) {
              FluxDir[Right][qp][var][idx(xdir, ydir)] =
                  Cons[Tidx(var, xdir, ydir)];
            }

#if NDIMS > 1
            if (ydir < yDim - 1) {
              FluxDir[Bottom][qp][var][xdir * yDim + ydir] =
                  Cons[Tidx(var, xdir, ydir)];
            }

            if (ydir > 0) {
              FluxDir[Top][qp][var][xdir * yDim + ydir] =
                  Cons[Tidx(var, xdir, ydir)];
            }
#endif
          }
        }
      }
    }
  }

  void WENO() {

    for (int qp = 0; qp < nqp; ++qp) {
      for (int var = 0; var < NumVar; ++var) {
        for (int xdir = XStart - 1; xdir < XEnd + 1; ++xdir) {
          for (int ydir = YStart; ydir < YEnd; ++ydir) {

            // Probably slow to do this pointer aliasing, but certainly makes
            // the code much easier to write
            double *Center = &Cons[Tidx(var, xdir, ydir)];

            double p1L = (-1.0 / 6.0) * (Center[-2]) +
                         (5.0 / 6.0) * (Center[-1]) + (1.0 / 3.0) * (Center[0]);

            double p1R = (1.0 / 3.0) * (Center[-2]) +
                         (-7.0 / 6.0) * (Center[-1]) +
                         (11.0 / 6.0) * (Center[0]);

            double p2L = (1.0 / 3.0) * (Center[-1]) +
                         (5.0 / 6.0) * (Center[0]) + (-1.0 / 6.0) * (Center[1]);

            double p2R = (-1.0 / 6.0) * (Center[-1]) +
                         (5.0 / 6.0) * (Center[0]) + (1.0 / 3.0) * (Center[1]);

            double p3L = (11.0 / 6.0) * (Center[0]) +
                         (-7.0 / 6.0) * (Center[1]) + (1.0 / 3.0) * (Center[2]);

            double p3R = (1.0 / 3.0) * (Center[0]) + (5.0 / 6.0) * (Center[1]) +
                         (-1.0 / 6.0) * (Center[2]);

            double Beta1 =
                (13.0 / 12.0) *
                    (std::pow(Center[-2] - 2.0 * Center[-1] + Center[0], 2)) +
                0.25 * (std::pow(
                           Center[-2] - 4.0 * Center[-1] + 3.0 * Center[0], 2));

            double Beta2 =
                (13.0 / 12.0) *
                    (std::pow(Center[-1] - 2.0 * Center[0] + Center[1], 2)) +
                0.25 * (std::pow(Center[-1] - Center[1], 2));

            double Beta3 =
                (13.0 / 12.0) *
                    (std::pow(Center[0] - 2.0 * Center[1] + Center[2], 2)) +
                0.25 * (std::pow(3.0 * Center[0] - 4.0 * Center[1] + Center[2],
                                 2));

            double eps = 1E-36;

            double w1L, w2L, w3L, w1R, w2R, w3R, wLSum, wRSum;

            w1L = 0.3 / (eps + Beta1);
            w1R = 0.1 / (eps + Beta1);

            w2L = 0.6 / (eps + Beta2);
            w2R = 0.6 / (eps + Beta2);

            w3L = 0.1 / (eps + Beta3);
            w3R = 0.3 / (eps + Beta3);

            wLSum = w1L + w2L + w3L;
            wRSum = w1R + w2R + w3R;

            w1L /= wLSum;
            w1R /= wRSum;

            w2L /= wLSum;
            w2R /= wRSum;

            w3L /= wLSum;
            w3R /= wRSum;

            // Populate the Fluxes, of each variable

            FluxDir[Left][qp][var][idx(xdir, ydir)] =
                w1L * p1L + w2L * p2L + w3L * p3L;
            FluxDir[Right][qp][var][idx(xdir, ydir)] =
                w1R * p1R + w2R * p2R + w3R * p3R;
          }
        }

#if NDIMS > 1
        if (ydir < yDim - 1) {
          FluxDir[Bottom][qp][var][xdir * yDim + ydir] =
              Cons[Tidx(var, xdir, ydir)];
        }

        if (ydir > 0) {
          FluxDir[Top][qp][var][xdir * yDim + ydir] =
              Cons[Tidx(var, xdir, ydir)];
        }
#endif
      }
    }
  }

  void HLL() {

    int i;
    int iPlus1;

    for (int Dim = 0; Dim < NDIMS * 2; Dim += 2) {
      for (int quadpoint = 0; quadpoint < nqp; ++quadpoint) {
        for (int xdir = 1; xdir < xDim - 1; ++xdir) {
          for (int ydir = 0; ydir < yDim; ++ydir) {

            i = idx(xdir, ydir);
            iPlus1 = (Dim == 0) ? idx(xdir + 1, ydir) : idx(xdir, ydir + 1);

            // FluxDir[0]

            double VelL = FluxDir[Dim + 1][quadpoint][Dim + 1][i] /
                          FluxDir[Dim + 1][quadpoint][Dens][i];
            // M_FluxR[i] /  D_FluxR[i];

            double VelR = FluxDir[Dim][quadpoint][Dim + 1][iPlus1] /
                          FluxDir[Dim][quadpoint][Dens][iPlus1];
            // double VelR = M_FluxL[i + 1] / D_FluxL[i + 1];

            // This Pressure assumes 1D, make sure to change later
            double PresL = (GAMMA - 1.0) *
                           (FluxDir[Dim + 1][quadpoint][Ener][i] -
                            VelL * FluxDir[Dim + 1][quadpoint][MomX][i] / 2.0);

            // double PresL = (gamma - 1.0) * (E_FluxR[i] - VelL * M_FluxR[i]
            // / 2.0);
            //
            double PresR = (GAMMA - 1.0) *
                           (FluxDir[Dim][quadpoint][Ener][iPlus1] -
                            VelR * FluxDir[Dim][quadpoint][MomX][iPlus1] / 2.0);

            // double PresR =
            //     (gamma - 1.0) * (E_FluxL[i + 1] - VelR * M_FluxL[i + 1]
            //     / 2.0);

            double CsL =
                std::sqrt(GAMMA * PresL / FluxDir[Dim + 1][quadpoint][Dens][i]);
            double CsR = std::sqrt(GAMMA * PresR /
                                   FluxDir[Dim][quadpoint][Dens][iPlus1]);
            // double CsL = std::sqrt(gamma * PresL / D_FluxR[i]);
            // double CsR = std::sqrt(gamma * PresR / D_FluxL[i + 1]);

            double SL = std::fmin(VelL - CsL, VelR - CsR);
            double SR = std::fmax(VelL + CsL, VelR + CsR);

            // HLL Wave switch, ugly but we both already know how it works so
            // ¯\_(ツ)_/¯
            if (0.0 <= SL) {
              Flux[quadpoint][Dens][i] = FluxDir[Dim + 1][quadpoint][MomX][i];
              // D_Flux[i] = M_FluxR[i];
              Flux[quadpoint][MomX][i] =
                  FluxDir[Dim + 1][quadpoint][MomX][i] * VelL + PresL;
              // M_Flux[i] = M_FluxR[i] * VelL + PresL;
              Flux[quadpoint][Ener][i] =
                  VelL * (FluxDir[Dim + 1][quadpoint][Ener][i] + PresL);
              // E_Flux[i] = VelL * (E_FluxR[i] + PresL);
            } else if (0.0 <= SR) {

              double LF1 = FluxDir[Dim + 1][quadpoint][MomX][i];
              // D_Flux[i] = M_FluxR[i];
              double LF2 = FluxDir[Dim + 1][quadpoint][MomX][i] * VelL + PresL;
              // M_Flux[i] = M_FluxR[i] * VelL + PresL;
              double LF3 =
                  VelL * (FluxDir[Dim + 1][quadpoint][Ener][i] + PresL);

              // double LF1 = M_FluxR[i];
              // double LF2 = M_FluxR[i] * VelL + PresL;
              // double LF3 = VelL * (E_FluxR[i] + PresL);

              double RF1 = FluxDir[Dim][quadpoint][MomX][iPlus1];
              // D_Flux[i] = M_FluxR[i];
              double RF2 = FluxDir[Dim][quadpoint][MomX][iPlus1] * VelL + PresL;
              // M_Flux[i] = M_FluxR[i] * VelL + PresL;
              double RF3 =
                  VelL * (FluxDir[Dim][quadpoint][Ener][iPlus1] + PresL);

              // double RF1 = M_FluxL[i + 1];
              // double RF2 = M_FluxL[i + 1] * VelR + PresR;
              // double RF3 = VelR * (E_FluxL[i + 1] + PresR);

              Flux[quadpoint][Dens][i] =
                  (SR * LF1 - SL * RF1 +
                   SL * SR *
                       (FluxDir[Dim][quadpoint][Dens][iPlus1] -
                        FluxDir[Dim + 1][quadpoint][Dens][i])) /
                  (SR - SL);

              Flux[quadpoint][MomX][i] =
                  (SR * LF2 - SL * RF2 +
                   SL * SR *
                       (FluxDir[Dim][quadpoint][MomX][iPlus1] -
                        FluxDir[Dim + 1][quadpoint][MomX][i])) /
                  (SR - SL);

              // M_Flux[i] =
              //     SR * LF2 - SL * RF2 + SL * SR * (M_FluxL[i + 1] -
              //     M_FluxR[i]);
              // M_Flux[i] /= (SR - SL);
              //
              Flux[quadpoint][Ener][i] =
                  (SR * LF3 - SL * RF3 +
                   SL * SR *
                       (FluxDir[Dim][quadpoint][Ener][iPlus1] -
                        FluxDir[Dim + 1][quadpoint][Ener][i])) /
                  (SR - SL);

              // E_Flux[i] =
              //     SR * LF3 - SL * RF3 + SL * SR * (E_FluxL[i + 1] -
              //     E_FluxR[i]);
              // E_Flux[i] /= (SR - SL);

              // D_Flux[i] =
              //     SR * LF1 - SL * RF1 + SL * SR * (D_FluxL[i + 1] -
              //     D_FluxR[i]);
              // D_Flux[i] /= (SR - SL);

              // M_Flux[i] =
              //     SR * LF2 - SL * RF2 + SL * SR * (M_FluxL[i + 1] -
              //     M_FluxR[i]);
              // M_Flux[i] /= (SR - SL);

              // E_Flux[i] =
              //     SR * LF3 - SL * RF3 + SL * SR * (E_FluxL[i + 1] -
              //     E_FluxR[i]);
              // E_Flux[i] /= (SR - SL);

            } else {

              Flux[quadpoint][Dens][i] = FluxDir[Dim][quadpoint][MomX][iPlus1];
              // D_Flux[i] = M_FluxR[i];
              Flux[quadpoint][MomX][i] =
                  FluxDir[Dim][quadpoint][MomX][iPlus1] * VelR + PresR;
              // M_Flux[i] = M_FluxR[i] * VelL + PresL;
              Flux[quadpoint][Ener][i] =
                  VelR * (FluxDir[Dim][quadpoint][Ener][iPlus1] + PresR);

              // D_Flux[i] = M_FluxL[i + 1];
              // M_Flux[i] = M_FluxL[i + 1] * VelR + PresR;
              // E_Flux[i] = VelR * (E_FluxL[i + 1] + PresR);
            }
          }
        }
      }
    }
  }

  void Recon() {
    for (int quad = 0; quad < nqp; ++quad) {
      for (int var = 0; var < NumVar; ++var) {
        for (int x = XStart; x < XEnd; ++x) {
          for (int y = YStart; y < YEnd; ++y) {
            Cons[Tidx(var, x, y)] -=
                *dt *
                (Flux[quad][var][idx(x, y)] - Flux[quad][var][idx(x - 1, y)]) /
                dx;
          }
        }
      }
    }
  }
};

#endif // FLUXCLASS_H_
