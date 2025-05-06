#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <cfenv>
#include <gsl/gsl_roots.h>
#include <iomanip>
#include <iostream>
#include "../include/QuadExactSolver.hpp"
#define TOLS (1.e-8)

using namespace std;
// typedef double realkind; 

const realkind Gamma = (5.0 / 3);
const realkind Sigma = (Gamma / (Gamma - 1.0));

realkind RelSpeed(realkind a, realkind b) { return (a - b) / (1 - a * b); }

struct IntParams {
  realkind S;
  realkind A;
};

struct RarefactionParams {
  struct IntParams Params;
  int sign;
  realkind v;
  realkind p;
};

struct ShockParams {
  int sign;
  realkind h, p, r, v, lor;
};

struct RR_Pstar_Params {
  realkind Lv, LA, LS, Lp;
  realkind Rv, RA, RS, Rp;
  realkind v0;
};

struct RS_Pstar_Params {
  realkind v0;
  realkind Lv, LA, LS, Lp;
  realkind Rh, Rp, Rr, Rlor, Rv;
};
struct SS_Pstar_Params {
  realkind Lh, Lp, Lr, Lv, Llor;
  realkind Rh, Rp, Rr, Rv, Rlor;
  realkind v0;
};
struct RareFactionSample_Params {
  realkind xi, pR, pL, SL, AL, vL, sign;
};

extern "C" {

static double Integral1(double P, void *pram) {
  realkind p = (realkind)P;
  struct IntParams *params = (struct IntParams *)pram;
  realkind S = params->S;
  realkind A = params->A;
  realkind rho = pow((p / S), 1.0 / Gamma);
  realkind H = 1.0 + Sigma * p / rho;
  realkind Cs = sqrt(Gamma * p / (H * rho));
  return (double)sqrt(H * H + A * A * (1 - Cs * Cs)) / ((H * H + A * A) * rho * Cs);
};

realkind RarefactionVx(realkind p, RarefactionParams *params) {

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
  double result, error;
  gsl_function Int;

  // Set up the integral
  Int.function = &Integral1;
  Int.params = &params->Params;

  gsl_integration_qag(&Int, params->p, p, 0, TOLS, 1000, 6, w, &result,
                      &error);

  gsl_integration_workspace_free(w);

  realkind B1 = .5 * log((1.0 + params->v) / (1.0 - params->v));
  return tanh(B1 + (realkind)params->sign * result);
};

realkind ux(realkind xi, realkind S, realkind press, realkind A, realkind sign) {
  realkind rho = pow(press / S, 1.0 / Gamma);
  realkind h = 1.0 + Sigma * press / rho;
  realkind cs = sqrt(Gamma * press / (h * rho));
  realkind a = cs * h;
  realkind b = sign * sqrt(A * A * (1.0 - cs * cs) + h * h);
  return (a - b * xi) / (a * xi - b);
};

double SampleRarefactionWave(double pressu, void *Parms) {
  realkind result;
  realkind pressure = (realkind ) pressu;
  struct RareFactionSample_Params *params =
      (struct RareFactionSample_Params *)Parms;
  struct IntParams IntPar = {params->SL, params->AL};
  struct RarefactionParams RarePar = {IntPar, (int)params->sign, params->vL,
                                      params->pL};

  result = ux(params->xi, params->SL, pressure, params->AL, params->sign) -
         RarefactionVx(pressure, &RarePar);

  return (double) result;
};

realkind Taub(realkind hA, realkind rho, realkind p, realkind pres) {
  realkind c_2 = (1.0 + (p - pres) / (pres * Sigma));
  realkind c_1 = -(p - pres) / (pres * Sigma);
  realkind c_0 = hA * (p - pres) / rho - hA * hA;
  realkind val, error = 1.0;

  if (fabs(c_2) < 1.0e-13) {
    val = -c_0 / c_1;
  } else {
    val = (-c_1 + sqrt(c_1 * c_1 - 4.0 * c_2 * c_0)) / (2.0 * c_2);
  }
  if (fabs(val - hA) < 1.e-15) {
    error = error / 0.0;
  }
  return val;
};

realkind J_sqr(realkind pres1, realkind pres2, realkind hA, realkind hB) {
  realkind val;
  if (fabs(hA - hB) < 1.0e-10) {
    val = Sigma * pres1 * pres2 / (hA * (hA - 1.0));
  } else {
    val = -Sigma * (pres1 - pres2) /
          (hA * (hA - 1.) / pres1 - hB * (hB - 1.0) / pres2);
  }

  if (val == val + 1.0) {
    cout << "landed here" << endl;
    cout << "WaveCrash" << endl;
    realkind error = 1 / val;
    error = 1.0 / val;
    error = error / 0.0;
  }
  return val;
};

realkind ShockSpeed(realkind lor, realkind r, realkind v, realkind J, int sign) {
  realkind D = r * lor;
  return (D * D * v + (realkind)sign * J * sqrt(J * J + D * D * (1.0 - v * v))) /
         (D * D + J * J);
}

realkind ShockVx(realkind p, ShockParams *params) {
  realkind hA = params->h;
  realkind hB = Taub(hA, params->r, params->p, p);
  realkind J2 = J_sqr(params->p, p, hA, hB);
  realkind J = sqrt(fabs(J2));
  realkind Vs = ShockSpeed(params->lor, params->r, params->v, J, params->sign);
  realkind Ws = 1.0 / sqrt(1.0 - Vs * Vs);

  return (hA * params->lor * params->v +
          (realkind)params->sign * Ws * (p - params->p) / J) /
         (hA * params->lor +
          (p - params->p) * ((realkind)params->sign * Ws * params->v / J +
                             1.0 / (params->r * params->lor)));
};

double DoubleRarefactionPstar( double pre, void *param) {
  realkind p = (realkind ) pre;
  realkind result;
  struct RR_Pstar_Params *params = (struct RR_Pstar_Params *)param;
  struct IntParams IntPar = {params->LS, params->LA};
  struct RarefactionParams Rare_Vx_Param = {IntPar, -1, params->Lv, params->Lp};

  realkind ux3 = RarefactionVx(p, &Rare_Vx_Param);

  IntPar.A = params->RA;
  IntPar.S = params->RS;
  Rare_Vx_Param.Params = IntPar;
  Rare_Vx_Param.sign = 1;
  Rare_Vx_Param.v = params->Rv;
  Rare_Vx_Param.p = params->Rp;

  realkind ux4 = RarefactionVx(p, &Rare_Vx_Param);

  realkind v13 = RelSpeed(params->Lv, ux3);
  realkind v24 = RelSpeed(params->Rv, ux4);

  result =  RelSpeed(v13, v24) - params->v0;
  return (double) result;
};

 double RareShockPstar(double pre, void *param) {
  realkind p = (realkind ) pre;
  realkind result;
  struct RS_Pstar_Params *params = (struct RS_Pstar_Params *)param;
  struct IntParams IntPar = {params->LS, params->LA};
  struct RarefactionParams Rare_Vx_Param = {IntPar, -1, params->Lv, params->Lp};

  realkind ux3 = RarefactionVx(p, &Rare_Vx_Param);
  struct ShockParams Shock_Pars;
  Shock_Pars.sign = 1;
  Shock_Pars.h = params->Rh;
  Shock_Pars.p = params->Rp;
  Shock_Pars.r = params->Rr;
  Shock_Pars.v = params->Rv;
  Shock_Pars.lor = params->Rlor;

  realkind ux4 = ShockVx(p, &Shock_Pars);

  realkind v13 = RelSpeed(params->Lv, ux3);
  realkind v64 = RelSpeed(params->Rv, ux4);

  result = RelSpeed(v13, v64) - params->v0;
  return (double) result;
};

double ShockShockPstar(double pre, void *param) {
  realkind p = (realkind ) pre;
  realkind result;
  struct SS_Pstar_Params *params = (struct SS_Pstar_Params *)param;

  struct ShockParams LShock = {
      -1, params->Lh, params->Lp, params->Lr, params->Lv, params->Llor,
  };
  struct ShockParams RShock = {
      1, params->Rh, params->Rp, params->Rr, params->Rv, params->Rlor,
  };

  realkind ux3 = ShockVx(p, &LShock);
  realkind ux4 = ShockVx(p, &RShock);

  realkind v13 = RelSpeed(params->Lv, ux3);
  realkind v64 = RelSpeed(params->Rv, ux4);
  result = RelSpeed(v13, v64) - params->v0;
  return (double) result;
};
}

class Wave {
public:
  realkind rho;
  realkind v;
  realkind vt;
  realkind p;
  string WaveType = "Default";
  realkind S;
  realkind A;
  realkind lor;
  realkind h;
  realkind ContactSpeed;
  realkind LeftB, RightB;
  void SetState(realkind state[4]) {
    rho = state[0];
    v = state[1];
    vt = state[2];
    p = state[3];

    SetEntropy();
    Setlor();
    SetHFromState();
    SetA();
  };

  void SampleRare(realkind xi, Wave *State, realkind sign, realkind Sample[4]) {
    realkind pmin = fmin(State->p, p);
    realkind pmax = fmax(State->p, p);
    realkind p_star;
    struct RareFactionSample_Params RareP = {
        xi, p, State->p, State->S, State->A, State->v, sign};
    int status;
    const int max_iter = 150;
    int iter = 0;

    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    gsl_function F;

    F.function = &SampleRarefactionWave;
    F.params = &RareP;
    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc(T);
    gsl_root_fsolver_set(s, &F, pmin, pmax);
    if(SampleRarefactionWave(pmin,&RareP)*SampleRarefactionWave(pmax, &RareP) >= 0.0){
      cout << "Failure Line 286" << endl;
      exit(0);
    }

    do {
      iter++;
      status = gsl_root_fsolver_iterate(s);
      p_star = gsl_root_fsolver_root(s);
      pmin = gsl_root_fsolver_x_lower(s);
      pmax = gsl_root_fsolver_x_upper(s);
      status = gsl_root_test_interval(pmin, pmax, 0, TOLS);

    } while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free(s);
    if (iter == max_iter) {
      cout << "Failed to converge before max iteration" << endl;
      exit(0);
    }

    Sample[3] = p_star;
    Sample[0] = pow(p_star / State->S, 1.0 / Gamma);
    realkind H = 1.0 + Sigma * p_star / Sample[0];
    // realkind ux(realkind xi, realkind S, realkind press, realkind A, realkind sign) {
    Sample[1] = ux(xi, State->S, p_star, State->A, sign);
    Sample[2] = State->A * sqrt((1.0 - Sample[1] * Sample[1]) /
                                (H * H + State->A * State->A));
  };

  void SetEntropy() { S = p / pow(rho, Gamma); };
  void Setlor() { lor = 1.0 / (sqrt(1.0 - vt * vt - v * v)); };
  void SetHFromState() { h = 1.0 + Sigma * p / rho; };

  void SetA() { A = lor * vt * h; };

  void FillWave(Wave *Neighbor) {
    if (WaveType == "Rarefaction") {
      rho = pow((p / Neighbor->S), 1.0 / Gamma);
      SetHFromState();
      SetEntropy();
      vt = sqrt(Neighbor->A * Neighbor->A * (1.0 - v * v) /
                (h * h + Neighbor->A * Neighbor->A));
      Setlor();
    } else if (WaveType == "Shock") {
      realkind hB = Taub(Neighbor->h, Neighbor->rho, Neighbor->p, p);
      realkind A = Neighbor->A;
      vt = A * sqrt((1.0 - v * v) / (hB * hB + A * A));
      rho = Sigma / (hB - 1.0) * p;
    } else {
      cout << "Unknown wavetype. \nExiting" << endl;
      exit(0);
    }
    SetEntropy();
    Setlor();
    SetHFromState();
  };

  void PrintWave() {
    cout << "[rho = " << rho << "] [V = " << v << "] [Vt = " << vt
         << "] [p = " << p << "]" << endl;
  };

  realkind SVel(Wave *state, realkind sS, realkind sign) {
    realkind SP_sqr = state->v * state->v + state->vt * state->vt;
    realkind r = pow(state->p / sS, 1.0 / Gamma);
    realkind h = 1.0 + Sigma * state->p / r;
    realkind cs = sqrt(Gamma * state->p / (h * r));
    realkind cs_sqr = cs * cs;
    return (state->v * (1.0 - cs_sqr) +
            sign * cs *
                sqrt((1.0 - SP_sqr) * (1.0 - SP_sqr * cs_sqr -
                                       state->v * state->v * (1.0 - cs_sqr)))) /
           (1.0 - SP_sqr * cs_sqr);
  };

  void Boundaries(Wave *StateL, Wave *StateR, realkind sign) {
    ContactSpeed = v;

    if (WaveType == "Rarefaction") {
      realkind s = StateL->p / (pow(StateL->rho, Gamma));
      realkind s2 = StateR->p / (pow(StateR->rho, Gamma));
      if (fabs(s - s2) > 1.e-10 * s) {
        cout << "Fail line 256";
        exit(0);
      }

      LeftB = SVel(StateL, s, sign);
      RightB = SVel(StateR, s, sign);

    } else if (WaveType == "Shock") {

      realkind hA = StateL->h;
      realkind pres = StateR->p;
      realkind hB = Taub(hA, StateL->rho, StateL->p, pres);
      realkind J2 = J_sqr(StateL->p, StateR->p, hA, hB);
      realkind J = sqrt(fabs(J2));
      realkind Vs = ShockSpeed(StateL->lor, StateL->rho, StateL->v, J, (int)sign);

      LeftB = Vs;
      RightB = Vs;

    } else {
      cout << "I don't know what I am :'(" << endl;
      exit(0);
    }
  };
};

class RiemannFan {
public:
  Wave WaveL;
  Wave WaveR;
  Wave Wave3;
  Wave Wave4;
  bool Reversed;

  void LoadStates(realkind StateL[4], realkind StateR[4]) {
    if (StateR[3] > StateL[3]) {
      Reversed = true;
      WaveR.SetState(StateL);
      WaveL.SetState(StateR);
    } else {
      Reversed = false;
      WaveL.SetState(StateL);
      WaveR.SetState(StateR);
    }
  };
  realkind FindD(const Wave &Left, const Wave &Right) const {
    realkind D;
    D = 4.0 * Gamma * Left.p;
    D *= ((Gamma - 1.0) * Right.p + Left.p) /
         pow((Gamma - 1.0) * (Left.p - Right.p), 2);
    D *= Right.h * (Right.p - Left.p) / Right.rho - Right.h * Right.h;
    return 1.0 - D;
  };

  realkind VaccuumLimit() const {
    realkind v1_x, v2_x;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    double result, error;
    gsl_function Int;

    // Set up the first integral
    struct IntParams Params = {WaveL.S, WaveL.A};
    Int.function = &Integral1;
    Int.params = &Params;

    /*
    gsl_integration_qags is known to work for singularities
    consider replacing the below if encounter difficulties
    */

    //  Assumed no singularities 61pt Gauss-Kronrod
    gsl_integration_qag(&Int, WaveL.p, 0, 0, TOLS, 1000, 6, w, &result,
                        &error);
    // gsl_integration_qags(&Int, WaveL.p, 0, 0, 1.0e-13, 1000, w, &result,
    //                      &error);
    v1_x = (realkind) tanh(result);

    // Set up the second integral
    Params.A = WaveR.A;
    Params.S = WaveR.S;

    gsl_integration_qag(&Int, 0, WaveR.p, 0, TOLS, 1000, 6, w, &result,
                        &error);
    // gsl_integration_qags(&Int, 0, WaveR.p, 0, 1.0e-13, 1000, w, &result,
    //                      &error);

    v2_x = (realkind) tanh(result);

    gsl_integration_workspace_free(w);

    return RelSpeed(v1_x, v2_x);
  };

  realkind DoubleRarefactionLimit() const {
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    double result, error;
    gsl_function Int;

    // Set up the integral
    struct IntParams Params = {WaveL.S, WaveL.A};
    Int.function = &Integral1;
    Int.params = &Params;

    gsl_integration_qag(&Int, WaveL.p, WaveR.p, 0, TOLS, 1000, 6, w, &result,
                        &error);

    gsl_integration_workspace_free(w);

    return (realkind) tanh(result);
  };

  realkind RareShockLimit() const {
    realkind D = FindD(WaveL, WaveR);
    realkind h3;
    realkind J2;
    realkind Vs;
    realkind Limit;

    h3 = (sqrt(D) - 1.0) * (Gamma - 1.0) * (WaveL.p - WaveR.p);
    h3 /= 2.0 * ((Gamma - 1.0) * WaveR.p + WaveL.p);

    J2 = -Sigma * (WaveL.p - WaveR.p);
    J2 /= h3 * (h3 - 1.0) / WaveL.p - WaveR.h * (WaveR.h - 1.0) / WaveR.p;

    Vs = pow(WaveR.rho * WaveR.lor, 2) * WaveR.v;
    Vs += sqrt(fabs(J2) * (J2 + pow(WaveR.rho * WaveR.lor, 2) *
                                    (1.0 - WaveR.v * WaveR.v)));
    Vs /= pow(WaveR.rho * WaveR.lor, 2) + J2;

    Limit = (WaveL.p - WaveR.p) * (1.0 - WaveR.v * Vs);
    Limit /= (Vs - WaveR.v) * (WaveR.h * WaveR.rho * WaveR.lor * WaveR.lor *
                                   (1.0 - WaveR.v * WaveR.v) +
                               WaveL.p - WaveR.p);

    return Limit;
  };

  void DoubleRarefactionStarValues(realkind v0) {
    realkind eps = 1.0e-10;
    realkind p_min = (WaveR.p + eps) * eps;
    realkind p_max = 2*WaveL.p;
    realkind v_star, p_star;
    int status;

    // cout << "V12_0 Values = " << v0 << endl;
    // if (p_min > p_max) {
    //   cout << "Waves facing the wrong way, 2 Rarefaction case" << endl;
    //   cout << "Terminating" << endl;
    //   exit(0);
    // }

    const int max_iter = 150;
    int iter = 0;

    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    gsl_function F;
    struct RR_Pstar_Params Parameters;
    Parameters.LA = WaveL.A;
    Parameters.LS = WaveL.S;
    Parameters.Lv = WaveL.v;
    Parameters.Lp = WaveL.p;

    Parameters.RA = WaveR.A;
    Parameters.RS = WaveR.S;
    Parameters.Rv = WaveR.v;
    Parameters.Rp = WaveR.p;

    Parameters.v0 = v0;

    F.function = &DoubleRarefactionPstar;
    F.params = &Parameters;

    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc(T);
    gsl_root_fsolver_set(s, &F, p_min, p_max);

    if(DoubleRarefactionPstar(p_min,&Parameters)*DoubleRarefactionPstar(p_max, &Parameters) >= 0.0){
      cout << "Failure Line 545" << endl;
      exit(0);
    }
    do {
      iter++;
      status = gsl_root_fsolver_iterate(s);
      p_star = gsl_root_fsolver_root(s);
      p_min = gsl_root_fsolver_x_lower(s);
      p_max = gsl_root_fsolver_x_upper(s);
      status = gsl_root_test_interval(p_min, p_max, 0, TOLS);

    } while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free(s);
    if (iter == max_iter) {
      cout << "Failed to converge before max iteration" << endl;
      exit(0);
    }

    struct IntParams Int = {WaveR.S, WaveR.A};
    struct RarefactionParams RareP = {Int, 1, WaveR.v, WaveR.p};

    v_star = RarefactionVx(p_star, &RareP);

    // cout << "P_star Value is: " << p_star << endl;
    // cout << "Converged in: " << iter << " iterations" << endl;
    // cout << "V_star Value is: " << v_star << endl;

    Wave3.v = v_star;
    Wave4.v = v_star;

    Wave3.p = p_star;
    Wave4.p = p_star;
  };

  void RareShockStarValues(realkind v12_0) {
    realkind eps = 1.0e-7;
    realkind p_min = WaveR.p - eps;
    realkind p_max = WaveL.p + eps;
    realkind v_star, p_star;
    int status;
    if (p_min > p_max) {
      cout << "Waves facing the wrong way, 2 Rarefaction case" << endl;
      cout << "Terminating" << endl;
      exit(0);
    }
    const int max_iter = 150;
    int iter = 0;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    gsl_function F;
    struct RS_Pstar_Params Parameters;

    Parameters.v0 = v12_0;

    Parameters.Lv = WaveL.v;
    Parameters.LA = WaveL.A;
    Parameters.LS = WaveL.S;
    Parameters.Lp = WaveL.p;

    Parameters.Rh = WaveR.h;
    Parameters.Rp = WaveR.p;
    Parameters.Rr = WaveR.rho;
    Parameters.Rlor = WaveR.lor;
    Parameters.Rv = WaveR.v;

    F.function = &RareShockPstar;
    F.params = &Parameters;

    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc(T);
    gsl_root_fsolver_set(s, &F, p_min, p_max);

    while(RareShockPstar(p_min,&Parameters)*RareShockPstar(p_max, &Parameters) >= 0.0){
      cout << "Failure Line 620" << endl;
      // exit(0);
      p_min *= .10;
      p_max *= 10.;
    }

    do {
      iter++;
      status = gsl_root_fsolver_iterate(s);
      p_star = gsl_root_fsolver_root(s);
      p_min = gsl_root_fsolver_x_lower(s);
      p_max = gsl_root_fsolver_x_upper(s);
      status = gsl_root_test_interval(p_min, p_max, 0.0, TOLS);

    } while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free(s);
    if (iter == max_iter) {
      cout << "Failed to converge before max iteration" << endl;
      exit(0);
    }
    struct ShockParams Shock_Par;
    Shock_Par.h = WaveR.h;
    Shock_Par.lor = WaveR.lor;
    Shock_Par.p = WaveR.p;
    Shock_Par.r = WaveR.rho;
    Shock_Par.sign = 1;
    Shock_Par.v = WaveR.v;
    v_star = ShockVx(p_star, &Shock_Par);

    Wave3.v = v_star;
    Wave4.v = v_star;

    Wave3.p = p_star;
    Wave4.p = p_star;

    // cout << "P_star Value is: " << p_star << endl;
    // cout << "Converged in: " << iter << " iterations" << endl;
    // cout << "V_star Value is: " << v_star << endl;
  };

  void TwoShockStarValues(realkind v12_0) {
    realkind v_star, p_star;
    realkind eps = 1.0e-10;
    // realkind p_min = WaveR.p, p_max = 10.0; // 1.0e11;
    realkind p_min = WaveL.p - eps;
    realkind p_max = 1.e12;
    int status;
    const int max_iter = 150;
    int iter = 0;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    gsl_function F;
    struct SS_Pstar_Params Parameters;
    Parameters.Lh = WaveL.h;
    Parameters.Lp = WaveL.p;
    Parameters.Lr = WaveL.rho;
    Parameters.Lv = WaveL.v;
    Parameters.Llor = WaveL.lor;

    Parameters.Rh = WaveR.h;
    Parameters.Rp = WaveR.p;
    Parameters.Rr = WaveR.rho;
    Parameters.Rv = WaveR.v;
    Parameters.Rlor = WaveR.lor;

    Parameters.v0 = v12_0;

    F.function = &ShockShockPstar;
    F.params = &Parameters;

    T = gsl_root_fsolver_brent;
    // T = gsl_root_fsolver_bisection;
    // T = gsl_root_fsolver_falsepos;
    s = gsl_root_fsolver_alloc(T);

    // cout << "Min Value = " << ShockShockPstar(p_min, &Parameters) << endl;
    // cout << "Max Value = " << ShockShockPstar(p_max, &Parameters) << endl;
    // exit(0);

    
    while(ShockShockPstar(p_min,&Parameters)*ShockShockPstar(p_max, &Parameters) >= 0.0){
      cout << "Failure Line 699" << endl;
      // exit(0);
      p_min *= .1;
    }
    gsl_root_fsolver_set(s, &F, p_min, p_max);
    do {
      iter++;
      status = gsl_root_fsolver_iterate(s);
      p_star = gsl_root_fsolver_root(s);
      p_min = gsl_root_fsolver_x_lower(s);
      p_max = gsl_root_fsolver_x_upper(s);
      status = gsl_root_test_interval(p_min, p_max, 0, TOLS);

    } while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free(s);
    if (iter == max_iter) {
      cout << "Failed to converge before max iteration" << endl;
      exit(0);
    }

    struct ShockParams Shock_Par;
    Shock_Par.h = WaveR.h;
    Shock_Par.lor = WaveR.lor;
    Shock_Par.p = WaveR.p;
    Shock_Par.r = WaveR.rho;
    Shock_Par.sign = 1;
    Shock_Par.v = WaveR.v;
    v_star = ShockVx(p_star, &Shock_Par);

    Wave3.v = v_star;
    Wave4.v = v_star;

    Wave3.p = p_star;
    Wave4.p = p_star;

    // cout << "P_star Value is: " << p_star << endl;
    // cout << "Converged in: " << iter << " iterations" << endl;
    // cout << "V_star Value is: " << v_star << endl;
  };

  void FindWaveTypes() {
    realkind v12_0 = RelSpeed(WaveL.v, WaveR.v);
    // cout << v12_0 << " The intermediate speed" << endl;
    // cout << "The realkind Rare Limit was " << DoubleRarefactionLimit() << endl;
    // cout << "The difference is " << DoubleRarefactionLimit() - v12_0 << endl;
    realkind Limit = 0.0; // VaccuumLimit();

    if (v12_0 <= Limit && false) {
      // This is an unneccessary check (hopefully). Remove when optimizing
      cout << "ERROR: Vaccuum case!" << endl;
      exit(0);
    } else {
      Limit = DoubleRarefactionLimit();
      if (v12_0 <= Limit || v12_0 - Limit < 1.e-15) {
        // This is the Two Rarefaction case
        Wave3.WaveType = "Rarefaction";
        Wave4.WaveType = "Rarefaction";
        DoubleRarefactionStarValues(v12_0);
      } else if (v12_0 <= RareShockLimit()) {
        // This is the One Rarefaction, One Shock case
        Wave3.WaveType = "Rarefaction";
        Wave4.WaveType = "Shock";
        RareShockStarValues(v12_0);
        // cout << "The limit was: " << RareShockLimit();
      } else {
        // This is the Two Shock case
        Wave3.WaveType = "Shock";
        Wave4.WaveType = "Shock";
        TwoShockStarValues(v12_0);
      }
    }
  };

  void CalculateIntermediateStates() {
    Wave3.FillWave(&WaveL);
    Wave4.FillWave(&WaveR);
  };

  void FanBoundaries() {
    Wave3.Boundaries(&WaveL, &Wave3, -1.0);
    Wave4.Boundaries(&Wave4, &WaveR, 1.0);
  };

  void SampleState(realkind xi, realkind State[4]) {
    if (Reversed) {
      xi = -xi;
    }
    if (xi < Wave3.LeftB) {
      // cout << "Left State" << endl;
      State[0] = WaveL.rho;
      State[1] = WaveL.v;
      State[2] = WaveL.vt;
      State[3] = WaveL.p;
    } else if (xi < Wave3.RightB) {
      // cout << "Rare State" << endl;
      Wave3.SampleRare(xi, &WaveL, -1.0, State);
      // cout << "In the Rarefaction Wave" << endl;
    } else if (xi < Wave3.ContactSpeed) {
      // cout << "Left Star" << endl;
      State[0] = Wave3.rho;
      State[1] = Wave3.v;
      State[2] = Wave3.vt;
      State[3] = Wave3.p;
    } else if (xi < Wave4.LeftB) {
      // cout << "Right Star" << endl;
      State[0] = Wave4.rho;
      State[1] = Wave4.v;
      State[2] = Wave4.vt;
      State[3] = Wave4.p;
    } else if (xi < Wave4.RightB) {
      // cout << "Right Rare" << endl;
      Wave4.SampleRare(xi, &WaveR, 1.0, State);
    } else {
      // cout << "Right State" << endl;
      State[0] = WaveR.rho;
      State[1] = WaveR.v;
      State[2] = WaveR.vt;
      State[3] = WaveR.p;
    }
  };
};

void SolveRiemannFlux(realkind StateL[4], realkind StateR[4], realkind Result[4], realkind T) {
  if (fabs(StateL[3] - StateR[3]) < 1.e-11){
    if (fabs(StateL[3] - StateR[3]) < 1.e-14){
      for (int i = 0; i < 4; ++i ){
        Result[i] = .5*(StateL[i] + StateR[i]);
      }
      // cout << "Stopped solve short" << endl;
      return;
    }
    double maxval = 0.0;
    for (int i = 0; i < 3; ++i){
      maxval = fmax( fabs(StateL[i] - StateR[i]), maxval);
    }
    if (maxval < 1.e-13){
      for (int i = 0; i < 4; ++i ){
        Result[i] = .5*(StateL[i] + StateR[i]);
      }
      // cout << "Stopped solve short" << endl;
      return;
    }
  }
  RiemannFan Problem;
  Problem.LoadStates(StateL, StateR);
  Problem.FindWaveTypes();
  Problem.CalculateIntermediateStates();
  Problem.FanBoundaries();
  Problem.SampleState(T, Result);
  // delete Problem;
}
void SolveShockTube(realkind StateL[4], realkind StateR[4], realkind Time) {
  RiemannFan Problem;
  realkind Result[4];
  realkind Dens[400], XVel[400], YVel[400], Pres[400];
  Problem.LoadStates(StateL, StateR);
  Problem.FindWaveTypes();
  Problem.CalculateIntermediateStates();
  Problem.FanBoundaries();
  for (int x = 0; x < 400; ++x) {
    Problem.SampleState((-.5 + (x / 399.0)) / Time, Result);
    Dens[x] = Result[0];
    XVel[x] = Result[1];
    YVel[x] = Result[2];
    Pres[x] = Result[3];
  }
  FILE *File1 = fopen("OutputData/Density.dat", "w");
  FILE *File2 = fopen("OutputData/VelocityX.dat", "w");
  FILE *File3 = fopen("OutputData/VelocityY.dat", "w");
  FILE *File4 = fopen("OutputData/Pressure.dat", "w");
  if (File1 && File2 && File3 && File4) {

    for (int i = 0; i < 400; i++) {

      fprintf(File1, "%.9Lg ", Dens[i]);
      fprintf(File2, "%.9Lg ", XVel[i]);
      fprintf(File3, "%.9Lg ", YVel[i]);
      fprintf(File4, "%.9Lg ", Pres[i]);
    }
    fprintf(File1, "\n");
    fprintf(File2, "\n");
    fprintf(File3, "\n");
    fprintf(File4, "\n");

    fclose(File1);
    fclose(File2);
    fclose(File3);
    fclose(File4);
  }
}

// int main(){
  // cout << setprecision(15);
// //   // RiemannFan Problem;

// //   // Test Case
// //   realkind StateL[4] = {1.0, .0, 0.99, 1000.0};
// //   realkind StateR[4] = {1.0, 0.0, 0.99, .01};
// //   realkind State[4];


// int main(){
//   feenableexcept(FE_INVALID);
// // realkind StateL[4] = {1.50617912044735, 0.469992817665391, 0.0, 0.417613895434584 };
// // realkind StateR[4] = {1.0, 0.0, 0.0, 0.00999999999999979 };
// realkind StateL[4] = {1.11471317639988, 0.507473466608499, 0, 0.507428210356333}; 
// realkind StateR[4] = {1, 0, 0, 0.00999999999999979}; 
// realkind Out[4];
// SolveRiemannFlux(StateL,  StateR,Out, -1.195631);
// for(int i = 0; i < 4; ++i){
//   cout << Out[i] << " " << endl;
// }
// return 0;
// }

//   // SolveShockTube(StateL, StateR, .4);
//   // Problem.LoadStates(StateL, StateR);
//   // Problem.FindWaveTypes();
//   // Problem.CalculateIntermediateStates();
//   // Problem.FanBoundaries();
//   // // Problem.SampleState(0.0, State);
//   // // cout << State[3] <<endl;
  
//   // cout << Problem.Wave4.rho << endl;
//   // cout << Problem.Wave4.v << endl;
//   // cout << Problem.Wave4.vt << endl;
//   // cout << Problem.Wave4.p << endl;


//   // void * blank;
//   // realkind test = SimpsonsRule(quad,0.5, 4.0, 1.e-12 , blank);
//   // cout << test << endl;
//   // test = NewtonIteration(quad, 0.10, 1.e-12 , blank);
//   // cout << test << endl;
//   return 0;


// }
