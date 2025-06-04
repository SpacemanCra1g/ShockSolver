#include <cmath>
#include <iomanip>
#include <iostream>
#define GAMMA (5.0 / 3.0)
using namespace std;

double SIGN(double x) { return (x >= 0.0) ? 1.0 : -1.0; }

class WaveFan {
public:
  double L_pRho, L_v, L_vt, L_p;
  double R_pRho, R_v, R_vt, R_p;
  double L_h;

  double L_cRho, L_m, L_mt, L_E;
  double R_cRho, R_m, R_mt, R_E;
  double R_h;
  double L_Cs, R_Cs;
  double Lam_RL, Lam_RR;
  double Lam_CL, Lam_CR;
  double SL, SR;
  double lamStar;
  WaveFan(double StateL[4], double StateR[4]) {
    // Fill Prim vars
    L_pRho = StateL[0];
    L_v = StateL[1];
    L_vt = StateL[2];
    L_p = StateL[3];

    R_pRho = StateR[0];
    R_v = StateR[1];
    R_vt = StateR[2];
    R_p = StateR[3];

    double v2, Lor, alpha;

    // Fill Con Vars
    v2 = pow(L_v, 2) + pow(L_vt, 2);

    Lor = 1.0 / sqrt(1.0 - v2);

    L_h = 1.0 + (GAMMA / (GAMMA - 1.0)) * L_p / L_pRho;

    alpha = L_pRho * L_h * Lor * Lor;

    L_cRho = L_pRho * Lor;
    L_m = L_v * alpha;
    L_mt = L_vt * alpha;
    L_E = alpha - L_p;

    v2 = pow(R_v, 2) + pow(R_vt, 2);

    Lor = 1.0 / sqrt(1.0 - v2);

    R_h = 1.0 + (GAMMA / (GAMMA - 1.0)) * R_p / R_pRho;

    alpha = R_pRho * R_h * Lor * Lor;

    R_cRho = R_pRho * Lor;
    R_m = R_v * alpha;
    R_mt = R_vt * alpha;
    R_E = alpha - R_p;

    // Fill Sound speed

    L_Cs = GAMMA * L_p / (L_h * L_pRho);
    R_Cs = GAMMA * R_p / (R_h * R_pRho);

    // Signal Speed
    double lor = 1.0 / sqrt(1.0 - L_v * L_v - L_vt * L_vt);
    double sroot = L_Cs / (lor * lor * (1.0 - L_Cs));
    Lam_CR = (L_v + sqrt(sroot * (1.0 - L_v * L_v + sroot))) / (1.0 + sroot);
    Lam_CL = (L_v - sqrt(sroot * (1.0 - L_v * L_v + sroot))) / (1.0 + sroot);

    lor = 1.0 / sqrt(1.0 - R_v * R_v - R_vt * R_vt);
    sroot = R_Cs / (lor * lor * (1.0 - R_Cs));
    Lam_RR = (R_v + sqrt(sroot * (1.0 - R_v * R_v + sroot))) / (1.0 + sroot);
    Lam_RL = (R_v - sqrt(sroot * (1.0 - R_v * R_v + sroot))) / (1.0 + sroot);

    SL = fmin(Lam_CL, Lam_RL);
    SR = fmax(Lam_CR, Lam_RR);

    double AL, AR, BL, BR, a, b, c;

    AL = SL * L_E - L_m;
    AR = SR * R_E - R_m;

    BL = L_m * (SL - L_v) - L_p;

    BR = R_m * (SR - R_v) - R_p;

    a = AR * SL - AL * SR;
    b = AL + BL * SR - AR - BR * SL;
    c = BR - BL;

    double scrh = -0.5 * (b + SIGN(b) * sqrt(b * b - 4.0 * a * c));
    lamStar = c / scrh;
  }
};

void SolveShockTube(double StateL[4], double StateR[4], double Time) {

  WaveFan Fan(StateL, StateR);
  cout << "SL = " << Fan.SL * Time + .5 << endl;
  cout << "SR = " << Fan.SR * Time + .5 << endl;
  cout << "LamStar = " << Fan.lamStar * Time + .5 << endl;
  // double SL, SS, SR;
  // double Result[4];
  // double Dens[400], XVel[400], YVel[400], Pres[400];

  // for (int x = 0; x < 400; ++x) {

  //   Dens[x] = Result[0];
  //   XVel[x] = Result[1];
  //   YVel[x] = Result[2];
  //   Pres[x] = Result[3];
  // }
  // FILE *File1 = fopen("OutputData/Density.dat", "w");
  // FILE *File2 = fopen("OutputData/VelocityX.dat", "w");
  // FILE *File3 = fopen("OutputData/VelocityY.dat", "w");
  // FILE *File4 = fopen("OutputData/Pressure.dat", "w");
  // if (File1 && File2 && File3 && File4) {

  //   for (int i = 0; i < 400; i++) {

  //     fprintf(File1, "%.9g ", Dens[i]);
  //     fprintf(File2, "%.9g ", XVel[i]);
  //     fprintf(File3, "%.9g ", YVel[i]);
  //     fprintf(File4, "%.9g ", Pres[i]);
  //   }
  //   fprintf(File1, "\n");
  //   fprintf(File2, "\n");
  //   fprintf(File3, "\n");
  //   fprintf(File4, "\n");

  //   fclose(File1);
  //   fclose(File2);
  //   fclose(File3);
  //   fclose(File4);
  // }
}

int main() {
  cout << setprecision(15);

  // Test Case
  double StateL[4] = {1.0, .0, 0.9, 1000.0};
  double StateR[4] = {1.0, 0.0, 0.9, .01};

  // SR Case
  // double StateL[4] = {1.0, .5, 0.0, 1.0};
  // double StateR[4] = {.125, 0.0, 0.0, .1};

  // 2R Case
  // double StateL[4] = {1.0, 0.0, 0.999, 1};
  // double StateR[4] = {.125, 0.5, 0.0, .1};

  //  2S Case
  // double StateL[4] = {1.0, 0.5, 0.0, 1};
  // double StateR[4] = {.125, 0.0, 0.9, .1};

  // All normal Forward States Work. Let's try reversing the waves

  // double StateR[4] = {1.0, .0, 0.9, 1000.0};
  // double StateL[4] = {1.0, 0.0, 0.9, .01};

  // SR Case
  // double StateR[4] = {1.0, .5, 0.0, 1.0};
  // double StateL[4] = {.125, 0.0, 0.0, .1};

  // 2R Case
  // double StateR[4] = {1.0, 0.0, 0.999, 1};
  // double StateL[4] = {.125, 0.5, 0.0, .1};

  //  2S Case
  // double StateR[4] = {1.0, 0.5, 0.0, 1};
  // double StateL[4] = {.125, 0.0, 0.9, .1};

  SolveShockTube(StateL, StateR, .4);
  // double Result[4];
  // double StateL[4] = {0.237488597667484, 0.333823183668201,
  // 0.950230347920938,
  //                     7.91488442703943};
  // double StateR[4] = {0.263094468682757, 0.334768150497356,
  // 0.948934888668851,
  //                     7.54624591778466};
  // SolveRiemannFlux(StateL, StateR, Result);
  // cout << "Result is: " << Result[0] << " " << Result[1] << " " << Result[2]
  //      << " " << Result[3] << " " << endl;
  return 0;
}
