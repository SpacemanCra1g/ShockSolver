#include <cmath>
#include <gsl/gsl_integration.h>
#include <iomanip>
#include <iostream>
#include <string>
#define GAMMA (5. / 3.)
#define SIGMA (GAMMA / (GAMMA - 1.0))
using namespace std;

struct StateStruct {
  double v, vt, rho, p;
  double m, mt, rhoCon, E;
  double h, A, S, Cs2;
  double lor, sign;
};

struct IntParams {
  double S;
  double A;
};

extern "C" {

double Integral1(double p, void *pram) {
  struct IntParams *params = (struct IntParams *)pram;
  double S = params->S;
  double A = params->A;
  double rho = pow((p / S), 1.0 / GAMMA);
  double H = 1.0 + SIGMA * p / rho;
  double Cs = sqrt(GAMMA * p / (H * rho));
  return sqrt(H * H + A * A * (1 - Cs * Cs)) / ((H * H + A * A) * rho * Cs);
};

double IteratePressure(double Pb, void *Parm) {
  struct StateStruct *params = (struct StateStruct *)Parm;
  string WaveType = (Pb <= params->p) ? "Rarefaction" : "Shock";
  double Vb;

  if (WaveType == "Shock") {

    double hB, J2, Vx, Ws;
    // Calculate the enthalpy behind using Taub adiabat
    double a = 1.0 + (params->p - Pb) / (SIGMA * Pb);
    double b = -(params->p - Pb) / (SIGMA * Pb);
    double c =
        params->h * (params->p - Pb) / params->rho - params->h * params->h;

    if (fabs(a) < 1.e-13) {
      hB = -c / b;
    } else {
      hB = .5 * (-b + sqrt(b * b - 4.0 * a * c)) / a;
    }
    // cout << "hB = " << hB << endl;
    // cout << "a = " << a << endl;
    // cout << "b = " << b << endl;
    // cout << "c = " << c << endl;

    // Calculate J2
    J2 = -SIGMA * (params->p - Pb) /
         (params->h * (params->h - 1.0) / params->p - hB * (hB - 1.0) / Pb);

    // Calculate Vx
    Vx = params->rhoCon * params->rhoCon * params->v;
    Vx += params->sign *
          sqrt(fabs(J2) * (J2 + params->rhoCon * params->rhoCon *
                                    (1.0 - params->v * params->v)));
    Vx /= params->rhoCon * params->rhoCon + J2;

    Ws = 1.0 / sqrt(1.0 - Vx * Vx);
    Vb = params->h * params->lor * params->v +
         Ws * (Pb - params->p) / sqrt(fabs(J2));
    Vb /= params->h * params->lor +
          (Pb - params->p) *
              (Ws * params->v / sqrt(fabs(J2)) + 1.0 / (params->rhoCon));

  } else {
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    double result, error;
    gsl_function Int;
    struct IntParams IntPar = {params->S, params->A};

    Int.function = &Integral1;
    Int.params = &IntPar;
    gsl_integration_qag(&Int, params->p, Pb, 0, 1.0e-12, 1000, 6, w, &result,
                        &error);

    gsl_integration_workspace_free(w);

    Vb = tanh(.5 * log((1.0 + params->v) / (1.0 - params->v)) +
              params->sign * result);
  }
  return Vb;
};
}

class State {
public:
  double v, vt, rho, p;
  double m, mt, rhoCon, E;
  double h, A, S, Cs2;
  double lor;
  double sign;
  struct StateStruct Params;
  void init(double state[4], double dir = 0.0) {
    rho = state[0];
    v = state[1];
    vt = state[2];
    p = state[3];
    lor = sqrt(1.0 / (1.0 - v * v - vt * vt));
    h = 1.0 + SIGMA * p / rho;
    Cs2 = GAMMA * p / (h * rho);

    double alpha = rho * h * lor * lor;

    rhoCon = rho * lor;
    m = v * alpha;
    mt = vt * alpha;
    E = alpha - p;
    sign = dir;

    A = lor * vt * h;
    S = p / pow(rho, GAMMA);
    Params.v = v;
    Params.vt = vt;
    Params.rho = rho;
    Params.p = p;
    Params.m = m;
    Params.mt = mt;
    Params.rhoCon = rhoCon;
    Params.E = E;
    Params.h = h;
    Params.A = A;
    Params.S = S;
    Params.Cs2 = Cs2;
    Params.lor = lor;
    Params.sign = sign;
  }
};

class WaveFan {
public:
  State StateL, StateR;
  bool Reversed;
  double SR, SL;
  bool Solved;
  State Solution;
  WaveFan(double Left[4], double Right[4]) {
    if (Left[3] > Right[3]) {
      StateL.init(Left, -1.0);
      StateR.init(Right, 1.0);
      Reversed = false;
    } else {
      StateR.init(Left, 1.0);
      StateL.init(Right, -1.0);
      Reversed = true;
    }
    Solved = false;
  };
  void HLLC_Edges() {
    if (Solved) {
      return;
    }

    double sroot = StateL.Cs2 / (StateL.lor * StateL.lor * (1.0 - StateL.Cs2));
    double Lam_CR =
      (StateL.v + sqrt(sroot * (1.0 - StateL.v * StateL.v + sroot))) /
        (1.0 + sroot);
    double Lam_CL =
        (StateL.v - sqrt(sroot * (1.0 - StateL.v * StateL.v + sroot))) /
        (1.0 + sroot);

    sroot = StateR.Cs2 / (StateR.lor * StateR.lor * (1.0 - StateR.Cs2));
    double Lam_RR =
        (StateR.v + sqrt(sroot * (1.0 - StateR.v * StateR.v + sroot))) /
        (1.0 + sroot);
    double Lam_RL =
        (StateR.v - sqrt(sroot * (1.0 - StateR.v * StateR.v + sroot))) /
        (1.0 + sroot);

    SL = fmin(Lam_CL, Lam_RL);
    SR = fmax(Lam_CR, Lam_RR);

    if (SL > 0.0) {
      double value[4] = {StateL.rho, StateL.v, StateL.vt, StateL.p};
      Solution.init(value);
      Solved = true;
    } else if (SR < 0.0) {
      double value[4] = {StateR.rho, StateR.v, StateR.vt, StateR.p};
      Solution.init(value);
      Solved = true;
    }
    return;
  };

  void Find_Pstar() {
    struct StateStruct LParams = {
        StateL.v,  StateL.vt,     StateL.rho, StateL.p,   StateL.m,
        StateL.mt, StateL.rhoCon, StateL.E,   StateL.h,   StateL.A,
        StateL.S,  StateL.Cs2,    StateL.lor, StateL.sign};
    struct StateStruct RParams = {
        StateR.v,  StateR.vt,     StateR.rho, StateR.p,   StateR.m,
        StateR.mt, StateR.rhoCon, StateR.E,   StateR.h,   StateR.A,
        StateR.S,  StateR.Cs2,    StateR.lor, StateR.sign};

    double checkvalue = 1.e-3;
    FILE *File1 = fopen("Pres.dat", "w");
    while (checkvalue < 200.0) {
      double test = IteratePressure(checkvalue, &LParams) -
                    IteratePressure(checkvalue, &RParams);
      // cout << "test value = " << test << endl;
      fprintf(File1, "%.9g ", test);
      checkvalue += .001;
      cout << "Check value = " << checkvalue << endl;
    }
    fclose(File1);
  };
};

int main() {
  cout << setprecision(15);
  double Left[4] = {1.0, .0, 0.9, 1000.0};
  double Right[4] = {1.0, 0.0, 0.9, .01};
  WaveFan Domain(Left, Right);
  Domain.Find_Pstar();
  return 0;
}


