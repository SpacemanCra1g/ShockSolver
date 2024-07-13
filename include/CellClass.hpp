#ifndef CELLCLASS_H_
#define CELLCLASS_H_

#include "GP_Kernel.hpp"
#include <string>

class FluxWall {

public:
  double *E_Values;
  double *MX_Values;
  double *MY_Values;
  double *D_Values;
  double *Vars[4];
  void FluxWall_init(int nqp) {
    E_Values = new double[nqp];
    MX_Values = new double[nqp];
    MY_Values = new double[nqp];
    D_Values = new double[nqp];
    Vars[0] = D_Values;
    Vars[1] = MX_Values;
    Vars[2] = MY_Values;
    Vars[3] = E_Values;
  }
};

class Cell {
public:
  double *DENS, *PRES, *XVEL, *YVEL, *MOMX, *MOMY, *ENERGY, gamma;
  double *Cs, dx, dy, *dt;
  Cell *LCell, *RCell, *TCell, *BCell;
  int x, y, nqp;
  FluxWall Walls[4];
  GP_Kernel *GP;
  double *p_R1_Stencil[20]; // 20 for NumVariable * Stencil length
  double R1_Stencil[20];    // Access via VarType * 5 + Stencil_IDX
  double UpdatedDens;
  double UpdatedMomX;
  double UpdatedMomY;
  double UpdatedE;

  // Defined VarConvert.cpp
  void Cons2Prims();
  void Prims2Cons();
  double GetPres();

  void SetUpCell(int i, int j, double *D, double *P, double *XV, double *YV,
                 double *MX, double *MY, double *E, double *CS, double dX,
                 double DY, double *Dt, double gam, GP_Kernel *Sol, int ngps) {

    DENS = D;
    PRES = P;
    XVEL = XV;
    YVEL = YV;
    MOMX = MX;
    MOMY = MY;
    ENERGY = E;
    Cs = CS;
    dx = dX;
    dy = DY;
    dt = Dt;
    gamma = gam;
    GP = Sol;
    nqp = ngps;
    x = i;
    y = j;
  }

  void Assign(std::string Var, double val) {
    if (Var == "DENS") {
      *DENS = val;
    } else if (Var == "PRES") {
      *PRES = val;
    } else if (Var == "XVEL") {
      *XVEL = val;
    } else if (Var == "YVEL") {
      *YVEL = val;
    }
  }

  void Assign_Stencils() {

    for (int i = 0; i < 4; ++i) {
      Walls[i].FluxWall_init(nqp);
    }

    p_R1_Stencil[0] = DENS;
    p_R1_Stencil[1 * 5] = MOMX;
    p_R1_Stencil[2 * 5] = MOMY;
    p_R1_Stencil[3 * 5] = ENERGY;

    p_R1_Stencil[0 * 5 + 1] = LCell->DENS;
    p_R1_Stencil[1 * 5 + 1] = LCell->MOMX;
    p_R1_Stencil[2 * 5 + 1] = LCell->MOMY;
    p_R1_Stencil[3 * 5 + 1] = LCell->ENERGY;

    p_R1_Stencil[0 * 5 + 2] = TCell->DENS;
    p_R1_Stencil[1 * 5 + 2] = TCell->MOMX;
    p_R1_Stencil[2 * 5 + 2] = TCell->MOMY;
    p_R1_Stencil[3 * 5 + 2] = TCell->ENERGY;

    p_R1_Stencil[0 * 5 + 3] = RCell->DENS;
    p_R1_Stencil[1 * 5 + 3] = RCell->MOMX;
    p_R1_Stencil[2 * 5 + 3] = RCell->MOMY;
    p_R1_Stencil[3 * 5 + 3] = RCell->ENERGY;

    p_R1_Stencil[0 * 5 + 4] = BCell->DENS;
    p_R1_Stencil[1 * 5 + 4] = BCell->MOMX;
    p_R1_Stencil[2 * 5 + 4] = BCell->MOMY;
    p_R1_Stencil[3 * 5 + 4] = BCell->ENERGY;
  }

  void Update_Stencils() {
    for (int i = 0; i < 20; ++i) {
      R1_Stencil[i] = *p_R1_Stencil[i];
      // std::cout << R1_Stencil[i] << std::endl;
    }
    // exit(0);
  }

  void GP_Reconstruction() {

    // The structure Walls[x].Y_Values[Z] with x = which wall, Y = var, Z = quad
    //

    Update_Stencils();

    for (int i = 0; i < 4; ++i) {
      Walls[0].Vars[i][0] = cblas_ddot(5, &R1_Stencil[5 * i], 1, GP->R1_g1R, 1);
      Walls[0].Vars[i][1] = cblas_ddot(5, &R1_Stencil[5 * i], 1, GP->R1_g2R, 1);

      Walls[1].Vars[i][0] = cblas_ddot(5, &R1_Stencil[5 * i], 1, GP->R1_g1L, 1);
      Walls[1].Vars[i][1] = cblas_ddot(5, &R1_Stencil[5 * i], 1, GP->R1_g2L, 1);

      Walls[2].Vars[i][0] = cblas_ddot(5, &R1_Stencil[5 * i], 1, GP->R1_g1U, 1);
      Walls[2].Vars[i][1] = cblas_ddot(5, &R1_Stencil[5 * i], 1, GP->R1_g2U, 1);

      Walls[3].Vars[i][0] = cblas_ddot(5, &R1_Stencil[5 * i], 1, GP->R1_g1D, 1);
      Walls[3].Vars[i][1] = cblas_ddot(5, &R1_Stencil[5 * i], 1, GP->R1_g2D, 1);
    }
  }

  void FindFlux(double &Dens, double &MOMX, double &MOMY, double &E,
                double Pres, int Dir) {
    double TD = Dens;
    double TX = MOMX;
    double TY = MOMY;
    double TE = E;

    if (Dir == 0) {
      Dens = MOMX;
      MOMY = TX * TY / TD;
      MOMX = TX * TX / TD + Pres;
      E = (TX / TD) * (TE + Pres);
    } else {
      Dens = MOMY;
      MOMX = TX * TY / TD;
      MOMY = TY * TY / TD + Pres;
      E = (TY / TD) * (TE + Pres);
    }
  }

  void HLLFlux(double &DL, double &XL, double &YL, double &EL, double &DR,
               double &XR, double &YR, double &ER, double SL, double SR,
               double PresL, double PresR, int Dir) {

    double DLT = DL;
    double XLT = XL;
    double YLT = YL;
    double ELT = EL;
    double DRT = DR;
    double XRT = XR;
    double YRT = YR;
    double ERT = ER;

    FindFlux(DL, XL, YL, EL, PresL, Dir);
    FindFlux(DR, XR, YR, ER, PresR, Dir);

    DL = (SR * DL - SL * DR + SR * SL * (DRT - DLT)) / (SR - SL);
    XL = (SR * XL - SL * XR + SR * SL * (XRT - XLT)) / (SR - SL);
    YL = (SR * YL - SL * YR + SR * SL * (YRT - YLT)) / (SR - SL);
    EL = (SR * EL - SL * ER + SR * SL * (ERT - ELT)) / (SR - SL);
  }

  void RiemannSolver() {
    FluxWall *LeftWall;
    FluxWall *RightWall;

    double MOMXL;
    double MOMXR;

    double EL;
    double ER;

    double DENSL;
    double DENSR;

    double MOMYL;
    double MOMYR;

    double VelL;
    double VelR;

    double PresL;

    double PresR;

    double CsL;
    double CsR;

    double SL;
    double SR;

    for (int Dir = 0; Dir < 3; Dir += 2) {
      for (int Gqp = 0; Gqp < nqp; ++Gqp) {

        LeftWall = &Walls[Dir];
        if (Dir == 0) {
          RightWall = &RCell->Walls[Dir + 1];
        } else {
          RightWall = &TCell->Walls[Dir + 1];
        }

        MOMXL = LeftWall->MX_Values[Gqp];
        MOMXR = RightWall->MX_Values[Gqp];

        EL = LeftWall->E_Values[Gqp];
        ER = RightWall->E_Values[Gqp];

        DENSL = LeftWall->D_Values[Gqp];
        DENSR = RightWall->D_Values[Gqp];

        MOMYL = LeftWall->MY_Values[Gqp];
        MOMYR = RightWall->MY_Values[Gqp];

        if (Dir == 0) {
          VelL = MOMXL / DENSL;
          VelR = MOMXR / DENSR;
        } else {
          VelL = MOMYL / DENSL;
          VelR = MOMYR / DENSR;
        }

        PresL = (gamma - 1.0) *
                (EL - 0.5 * (std::pow(MOMXL, 2) + std::pow(MOMYL, 2)) / DENSL);

        PresR = (gamma - 1.0) *
                (ER - 0.5 * (std::pow(MOMXR, 2) + std::pow(MOMYR, 2)) / DENSR);

        // Pres = (gamma - 1.0) * (E - .5 * (MOMX ^ 2 + MOMY ^ 2) / DENS);

        CsL = std::sqrt(gamma * PresL / DENSL);
        CsR = std::sqrt(gamma * PresR / DENSR);

        SL = std::fmin(VelL - CsL, VelR - CsR);
        SR = std::fmax(VelL + CsL, VelR + CsR);

        // HLL Wave switch, ugly but we both already know how it works so
        // ¯\_(ツ)_/¯
        if (0.0 <= SL) {
          FindFlux(DENSL, MOMXL, MOMYL, EL, PresL, Dir);
          LeftWall->D_Values[Gqp] = DENSL;
          LeftWall->MX_Values[Gqp] = MOMXL;
          LeftWall->MY_Values[Gqp] = MOMYL;
          LeftWall->E_Values[Gqp] = EL;

          RightWall->D_Values[Gqp] = DENSL;
          RightWall->MX_Values[Gqp] = MOMXL;
          RightWall->MY_Values[Gqp] = MOMYL;
          RightWall->E_Values[Gqp] = EL;
        } else if (0.0 <= SR) {

          HLLFlux(DENSL, MOMXL, MOMYL, EL, DENSR, MOMXR, MOMYR, ER, SL, SR,
                  PresL, PresR, Dir);
          LeftWall->D_Values[Gqp] = DENSL;
          LeftWall->MX_Values[Gqp] = MOMXL;
          LeftWall->MY_Values[Gqp] = MOMYL;
          LeftWall->E_Values[Gqp] = EL;

          RightWall->D_Values[Gqp] = DENSL;
          RightWall->MX_Values[Gqp] = MOMXL;
          RightWall->MY_Values[Gqp] = MOMYL;
          RightWall->E_Values[Gqp] = EL;
        } else {
          FindFlux(DENSR, MOMXR, MOMYR, ER, PresR, Dir);
          LeftWall->D_Values[Gqp] = DENSR;
          LeftWall->MX_Values[Gqp] = MOMXR;
          LeftWall->MY_Values[Gqp] = MOMYR;
          LeftWall->E_Values[Gqp] = ER;

          RightWall->D_Values[Gqp] = DENSR;
          RightWall->MX_Values[Gqp] = MOMXR;
          RightWall->MY_Values[Gqp] = MOMYR;
          RightWall->E_Values[Gqp] = ER;
        }
      }
    }
  }

  void FluxRecon() {
    UpdatedDens = 0.0;
    UpdatedE = 0.0;
    UpdatedMomX = 0.0;
    UpdatedMomY = 0.0;
    for (int i = 0; i < 2; ++i) {
      UpdatedDens += .5 * (Walls[0].D_Values[i] - Walls[2].D_Values[i]) / dx +
                     .5 * (Walls[1].D_Values[i] - Walls[3].D_Values[i]) / dy;

      UpdatedE += .5 * (Walls[0].E_Values[i] - Walls[2].E_Values[i]) / dx +
                  .5 * (Walls[1].E_Values[i] - Walls[3].E_Values[i]) / dy;

      UpdatedMomX += .5 * (Walls[0].MX_Values[i] - Walls[2].MX_Values[i]) / dx +
                     .5 * (Walls[1].MX_Values[i] - Walls[3].MX_Values[i]) / dy;

      UpdatedMomY += .5 * (Walls[0].MY_Values[i] - Walls[2].MY_Values[i]) / dx +
                     .5 * (Walls[1].MY_Values[i] - Walls[3].MY_Values[i]) / dy;
    }

    UpdatedDens = *DENS - *dt * UpdatedDens;
    UpdatedMomX = *MOMX - *dt * UpdatedMomX;
    UpdatedMomY = *MOMY - *dt * UpdatedMomY;
    UpdatedE = *ENERGY - *dt * UpdatedE;
  }
  void AcceptUpdate() {
    *DENS = UpdatedDens;
    *MOMX = UpdatedMomX;
    *MOMY = UpdatedMomY;
    *ENERGY = UpdatedE;
  }
};

#endif // CELLCLASS_H_
