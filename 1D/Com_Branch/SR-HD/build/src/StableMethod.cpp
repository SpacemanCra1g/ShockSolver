#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
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
  double lor, sign, otherP;
};

struct IntParams {
  double S;
  double A;
};

struct PressureStarRoot_Params {
  struct StateStruct LParams, RParams;
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
  if (fabs(Pb - params->p) < 1.e-13) {
    return params->v;
  }

  string WaveType;
  if (Pb < fmin(params->p, params->otherP)) {
    WaveType = "Rarefaction";
  } else if (Pb > fmax(params->p, params->otherP)) {
    WaveType = "Shock";
  } else if (params->sign > 0.0) {
    WaveType = "Shock";
  } else {
    WaveType = "Rarefaction";
  }

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

    if (std::isinf(J2)) {
      J2 = 1.e11;
    }
    // cout << "J2 = " << J2 << endl;
    // cout << "Pressure diff = " << Pb - params->p << endl;
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
  } else if (WaveType == "Rarefaction") {
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    double result, error;
    gsl_function Int;
    struct IntParams IntPar = {params->S, params->A};

    Int.function = &Integral1;
    Int.params = &IntPar;
    gsl_integration_qag(&Int, params->p, Pb, 0, 1.0e-8, 1000, 1, w, &result,
                        &error);

    gsl_integration_workspace_free(w);

    Vb = tanh(.5 * log((1.0 + params->v) / (1.0 - params->v)) +
              params->sign * result);
  } else {
    cout << "Unknown wave\n" << "Exiting" << endl;
  }
  return Vb;
};

double PressureStarRoot(double pStar, void *params) {
  struct PressureStarRoot_Params *Par =
      (struct PressureStarRoot_Params *)params;
  return IteratePressure(pStar, &Par->LParams) -
         IteratePressure(pStar, &Par->RParams);
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
    // cout << "Checkpoint 0" << endl;
    rho = state[0];
    v = state[1];
    vt = state[2];
    p = state[3];
    // cout << "Checkpoint 0.1" << endl;
    cout << -v * v - vt * vt << endl;
    lor = sqrt(1.0 / (1.0 - v * v - vt * vt));
    // cout << "Checkpoint 0.2" << endl;
    h = 1.0 + SIGMA * p / rho;
    // cout << "Checkpoint 0.3" << endl;
    Cs2 = GAMMA * p / (h * rho);

    // cout << "Checkpoint 1" << endl;
    double alpha = rho * h * lor * lor;

    rhoCon = rho * lor;
    m = v * alpha;
    mt = vt * alpha;
    E = alpha - p;
    sign = dir;

    A = lor * vt * h;
    S = p / pow(rho, GAMMA);
    Params.v = v;
    // cout << "Checkpoint 2" << endl;
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
    // cout << "Checkpoint 3" << endl;
    Params.lor = lor;
    Params.sign = sign;
    // cout << "Checkpoint 4" << endl;
  }
};

class WaveFan {
public:
  State StateL, StateR;
  bool Reversed;
  double SR, SL;
  bool Solved;
  State Solution;
  double p_star;
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
  void HLLC_Edges(double Out[4]) {
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
      Out[0] = StateL.rho;
      Out[1] = StateL.v;
      Out[2] = StateL.vt;
      Out[3] = StateL.p;
      Solved = true;
    } else if (SR < 0.0) {
      Out[0] = StateR.rho;
      Out[1] = StateR.v;
      Out[2] = StateR.vt;
      Out[3] = StateR.p;
      Solved = true;
    }
    return;
  };

  void Find_Pstar() {
    double p_min, p_max;
    struct StateStruct LParams = {
        StateL.v,  StateL.vt,     StateL.rho, StateL.p,    StateL.m,
        StateL.mt, StateL.rhoCon, StateL.E,   StateL.h,    StateL.A,
        StateL.S,  StateL.Cs2,    StateL.lor, StateL.sign, StateR.p};
    struct StateStruct RParams = {
        StateR.v,  StateR.vt,     StateR.rho, StateR.p,    StateR.m,
        StateR.mt, StateR.rhoCon, StateR.E,   StateR.h,    StateR.A,
        StateR.S,  StateR.Cs2,    StateR.lor, StateR.sign, StateL.p};

    double checkvalue = 1.e-1;
    // FILE *File1 = fopen("Pres.dat", "w");
    // while (checkvalue < 200.0) {
    double test = IteratePressure(checkvalue, &LParams) -
                  IteratePressure(checkvalue, &RParams);

    // double test = IteratePressure(checkvalue, &RParams);

    // cout << "test value = " << test << endl;
    // exit(0);
    // fprintf(File1, "%.9g ", test);
    // checkvalue += .001;
    if (test > 0.0) {
      p_min = checkvalue;
      p_max = checkvalue;
      while (test > 0.0) {
        p_max *= 10.0;
        test =
            IteratePressure(p_max, &LParams) - IteratePressure(p_max, &RParams);
      }
    } else {
      p_min = checkvalue;
      p_max = checkvalue;
      while (test < 0.0) {
        p_min *= 0.1;
        test =
            IteratePressure(p_min, &LParams) - IteratePressure(p_min, &RParams);
      }
    }

    int status;
    const int max_iter = 150;
    int iter = 0;

    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    gsl_function F;
    struct PressureStarRoot_Params Params = {LParams, RParams};

    F.function = &PressureStarRoot;
    F.params = &Params;
    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc(T);
    gsl_root_fsolver_set(s, &F, p_min, p_max);

    test = IteratePressure(p_min, &LParams) - IteratePressure(p_min, &RParams);
    // cout << "p_min = " << p_min << " Value is = " << test << endl;

    test = IteratePressure(p_max, &LParams) - IteratePressure(p_max, &RParams);
    // cout << "p_max = " << p_max << " Value is = " << test << endl;

    do {
      iter++;
      status = gsl_root_fsolver_iterate(s);
      p_star = gsl_root_fsolver_root(s);
      p_min = gsl_root_fsolver_x_lower(s);
      p_max = gsl_root_fsolver_x_upper(s);
      status = gsl_root_test_interval(p_min, p_max, 0, 1.e-10);

    } while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free(s);
    if (iter == max_iter) {
      cout << "Failed to converge before max iteration" << endl;
      exit(0);
    }

    // cout << "p_star = " << p_star << endl;
    // cout << "Converged in " << iter << " iterations." << endl;
  };

  void HLLCMethod(double Out[4]) {
    struct StateStruct RParams = {
        StateR.v,  StateR.vt,     StateR.rho, StateR.p,    StateR.m,
        StateR.mt, StateR.rhoCon, StateR.E,   StateR.h,    StateR.A,
        StateR.S,  StateR.Cs2,    StateR.lor, StateR.sign, StateL.p};

    struct StateStruct LParams = {
        StateL.v,  StateL.vt,     StateL.rho, StateL.p,    StateL.m,
        StateL.mt, StateL.rhoCon, StateL.E,   StateL.h,    StateL.A,
        StateL.S,  StateL.Cs2,    StateL.lor, StateL.sign, StateR.p};

    double v_star = IteratePressure(p_star, &LParams);

    State *state;

    if (v_star < 0.0) {
      state = &StateR;
    } else {
      state = &StateL;
    }
    if (p_star > fmax(StateL.p, StateR.p) ||
        (p_star > fmin(StateL.p, StateR.p) && v_star < 0.0)) {
      // Double Shock wave or Rare-Shock and we're in U*right
      double hB;
      double a = 1.0 + (state->p - p_star) / (SIGMA * p_star);
      double b = -(state->p - p_star) / (SIGMA * p_star);
      double c =
          state->h * (state->p - p_star) / state->rho - state->h * state->h;

      if (fabs(a) < 1.e-13) {
        hB = -c / b;
      } else {
        hB = .5 * (-b + sqrt(b * b - 4.0 * a * c)) / a;
      }

      double vt_c = state->A * sqrt((1.0 - v_star * v_star) /
                                    (hB * hB + state->A * state->A));

      double rho_c = SIGMA / (hB - 1.0) * p_star;
      Out[0] = rho_c;
      Out[1] = v_star;
      Out[2] = vt_c;
      Out[3] = p_star;
    } else {
      double rho_c = pow((p_star / state->S), 1.0 / GAMMA);
      double h_c = 1.0 + SIGMA * p_star / rho_c;
      double vt_c = sqrt(state->A * state->A * (1.0 - v_star * v_star) /
                         (h_c * h_c + state->A * state->A));
      Out[0] = rho_c;
      Out[1] = v_star;
      Out[2] = vt_c;
      Out[3] = p_star;
    }

    if (isnan(Out[3])) {
      cout << "Pressure Nan";
      exit(1);
    } else if (Out[3] <= 0.0) {
      cout << "Pressure Negative";
      exit(1);
    }
  };
};

void ExactHLL(double Left[4], double Right[4], double Out[4]) {
  WaveFan Domain(Left, Right);
  Domain.HLLC_Edges(Out);
  if (Domain.Solved) {
    return;
  }
  Domain.Find_Pstar();
  Domain.HLLCMethod(Out);
};

// int main() {
//   cout << setprecision(15);
//   double Out[4];
//   // double Left[4] = {1.0, .0, 0.9, 1000.0};
//   // double Right[4] = {1.0, 0.0, 0.9, .01};

//   double Left[4] = {0.156441069754414, 0.325321489139702, 0.941749072682333,
//                     3.29117532785489};
//   double Right[4] = {0.278903639690685, 0.362517872988394, 0.939646219630043,
//                      3.23680637291663};

//   // double Left[4] = {1.0, .0, 0.0, 1.0};
//   // double Right[4] = {.125, 0.5, 0.0, .1};

//   // WaveFan Domain(Left, Right);
//   ExactHLL(Left, Right, Out);
//   cout << Out[0] << " " << Out[1] << " " << Out[2] << " " << Out[3] << endl;

//   return 0;
// }
