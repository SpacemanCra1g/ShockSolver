#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <iomanip>
#include <iostream>

using namespace std;

const double Gamma = 5.0 / 3;
const double Sigma = Gamma / (Gamma - 1.0);

double RelSpeed(double a, double b) { return (a - b) / (1 - a * b); }

struct IntParams {
  double S;
  double A;
};

struct RarefactionParams {
  struct IntParams Params;
  int sign;
  double v;
  double p;
};

struct ShockParams {
  int sign;
  double h, p, r, v, lor;
};

struct RR_Pstar_Params {
  double Lv, LA, LS, Lp;
  double Rv, RA, RS, Rp;
  double v0;
};

struct RS_Pstar_Params {
  double v0;
  double Lv, LA, LS, Lp;
  double Rh, Rp, Rr, Rlor, Rv;
};
struct SS_Pstar_Params {
  double Lh, Lp, Lr, Lv, Llor;
  double Rh, Rp, Rr, Rv, Rlor;
  double v0;
};

extern "C" {

double Integral1(double p, void *pram) {
  struct IntParams *params = (struct IntParams *)pram;
  double S = params->S;
  double A = params->A;
  double rho = pow((p / S), 1.0 / Gamma);
  double H = 1.0 + Sigma * p / rho;
  double Cs = sqrt(Gamma * p / (H * rho));
  return sqrt(H * H + A * A * (1 - Cs * Cs)) / ((H * H + A * A) * rho * Cs);
};

double RarefactionVx(double p, RarefactionParams *params) {

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
  double result, error;
  gsl_function Int;

  // Set up the integral
  Int.function = &Integral1;
  Int.params = &params->Params;

  gsl_integration_qag(&Int, params->p, p, 0, 1.0e-12, 1000, 6, w, &result,
                      &error);

  gsl_integration_workspace_free(w);

  double B1 = .5 * log((1.0 + params->v) / (1.0 - params->v));
  return tanh(B1 + params->sign * result);
};

double Taub(double hA, double rho, double p, double pres) {
  double c_2 = (1.0 + (p - pres) / (pres * Sigma));
  double c_1 = -(p - pres) / (pres * Sigma);
  double c_0 = hA * (p - pres) / rho - hA * hA;

  if (fabs(c_2) < 1.0e-13) {
    return -c_0 / c_1;
  } else {
    return (-c_1 + sqrt(c_1 * c_1 - 4.0 * c_2 * c_0)) / (2.0 * c_2);
  }
};

double J_sqr(double pres1, double pres2, double hA, double hB) {
  double val = -Sigma * (pres1 - pres2) /
               (hA * (hA - 1.) / pres1 - hB * (hB - 1.0) / pres2);
  return val;
};

double ShockSpeed(double lor, double r, double v, double J, int sign) {
  double D = r * lor;
  return (D * D * v + (double)sign * J * sqrt(J * J + D * D * (1.0 - v * v))) /
         (D * D + J * J);
}

double ShockVx(double p, ShockParams *params) {
  double hA = params->h;
  double hB = Taub(hA, params->r, params->p, p);
  double J2 = J_sqr(params->p, p, hA, hB);
  double J = sqrt(fabs(J2));
  double Vs = ShockSpeed(params->lor, params->r, params->v, J, params->sign);
  double Ws = 1.0 / sqrt(1.0 - Vs * Vs);

  return (hA * params->lor * params->v +
          (double)params->sign * Ws * (p - params->p) / J) /
         (hA * params->lor +
          (p - params->p) * ((double)params->sign * Ws * params->v / J +
                             1.0 / (params->r * params->lor)));
};

double DoubleRarefactionPstar(double p, void *param) {
  struct RR_Pstar_Params *params = (struct RR_Pstar_Params *)param;
  struct IntParams IntPar = {params->LS, params->LA};
  struct RarefactionParams Rare_Vx_Param = {IntPar, -1, params->Lv, params->Lp};

  double ux3 = RarefactionVx(p, &Rare_Vx_Param);

  IntPar.A = params->RA;
  IntPar.S = params->RS;
  Rare_Vx_Param.Params = IntPar;
  Rare_Vx_Param.sign = 1;
  Rare_Vx_Param.v = params->Rv;
  Rare_Vx_Param.p = params->Rp;

  double ux4 = RarefactionVx(p, &Rare_Vx_Param);

  double v13 = RelSpeed(params->Lv, ux3);
  double v24 = RelSpeed(params->Rv, ux4);

  return RelSpeed(v13, v24) - params->v0;
};

double RareShockPstar(double p, void *param) {
  struct RS_Pstar_Params *params = (struct RS_Pstar_Params *)param;
  struct IntParams IntPar = {params->LS, params->LA};
  struct RarefactionParams Rare_Vx_Param = {IntPar, -1, params->Lv, params->Lp};

  double ux3 = RarefactionVx(p, &Rare_Vx_Param);
  struct ShockParams Shock_Pars;
  Shock_Pars.sign = 1;
  Shock_Pars.h = params->Rh;
  Shock_Pars.p = params->Rp;
  Shock_Pars.r = params->Rr;
  Shock_Pars.v = params->Rv;
  Shock_Pars.lor = params->Rlor;

  double ux4 = ShockVx(p, &Shock_Pars);

  double v13 = RelSpeed(params->Lv, ux3);
  double v64 = RelSpeed(params->Rv, ux4);

  return RelSpeed(v13, v64) - params->v0;
};

double ShockShockPstar(double p, void *param) {
  struct SS_Pstar_Params *params = (struct SS_Pstar_Params *)param;

  struct ShockParams LShock = {
      -1, params->Lh, params->Lp, params->Lr, params->Lv, params->Llor,
  };
  struct ShockParams RShock = {
      1, params->Rh, params->Rp, params->Rr, params->Rv, params->Rlor,
  };

  double ux3 = ShockVx(p, &LShock);
  double ux4 = ShockVx(p, &RShock);

  double v13 = RelSpeed(params->Lv, ux3);
  double v64 = RelSpeed(params->Rv, ux4);
  return RelSpeed(v13, v64) - params->v0;
};
}

class Wave {
public:
  double rho;
  double v;
  double vt;
  double p;
  string WaveType = "Default";
  double S;
  double A;
  double lor;
  double h;
  double ContactSpeed;
  double LeftB, RightB;
  void SetState(double state[4]) {
    rho = state[0];
    v = state[1];
    vt = state[2];
    p = state[3];

    SetEntropy();
    Setlor();
    SetHFromState();
    SetA();
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
      double hB = Taub(Neighbor->h, Neighbor->rho, Neighbor->p, p);
      double A = Neighbor->A;
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

  double SVel(Wave *state, double sS, double sign) {
    double SP_sqr = state->v * state->v + state->vt * state->vt;
    double r = pow(state->p / sS, 1.0 / Gamma);
    double h = 1.0 + Sigma * state->p / r;
    double cs = sqrt(Gamma * state->p / (h * r));
    double cs_sqr = cs * cs;
    return (state->v * (1.0 - cs_sqr) +
            sign * cs *
                sqrt((1.0 - SP_sqr) * (1.0 - SP_sqr * cs_sqr -
                                       state->v * state->v * (1.0 - cs_sqr)))) /
           (1.0 - SP_sqr * cs_sqr);
  };

  void Boundaries(Wave *StateL, Wave *StateR, double sign) {
    ContactSpeed = v;

    if (WaveType == "Rarefaction") {
      double s = StateL->p / (pow(StateL->rho, Gamma));
      double s2 = StateR->p / (pow(StateR->rho, Gamma));
      if (fabs(s - s2) > 1.e-10 * s) {
        cout << "Fail line 256";
        exit(0);
      }

      LeftB = SVel(StateL, s, sign);
      RightB = SVel(StateR, s, sign);

    } else if (WaveType == "Shock") {

      double hA = StateL->h;
      double pres = StateR->p;
      double hB = Taub(hA, StateL->rho, StateL->p, pres);
      cout << "Test:: hB == " << hA << endl;
      double J2 = J_sqr(StateL->p, StateR->p, hA, hB);
      double J = sqrt(fabs(J2));
      double Vs = ShockSpeed(StateL->lor, StateL->rho, StateL->v, J, (int)sign);

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

  void LoadStates(double StateL[4], double StateR[4]) {
    WaveL.SetState(StateL);
    WaveR.SetState(StateR);
  };
  double FindD(const Wave &Left, const Wave &Right) const {
    double D;
    D = 4.0 * Gamma * Left.p;
    D *= ((Gamma - 1.0) * Right.p + Left.p) /
         pow((Gamma - 1.0) * (Left.p - Right.p), 2);
    D *= Right.h * (Right.p - Left.p) / Right.rho - Right.h * Right.h;
    return 1 - D;
  };

  double VaccuumLimit() const {
    double v1_x, v2_x;
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
    gsl_integration_qag(&Int, WaveL.p, 0, 0, 1.0e-12, 1000, 6, w, &result,
                        &error);
    v1_x = tanh(result);

    // Set up the second integral
    Params.A = WaveR.A;
    Params.S = WaveR.S;

    gsl_integration_qag(&Int, 0, WaveR.p, 0, 1.0e-12, 1000, 6, w, &result,
                        &error);

    v2_x = tanh(result);

    gsl_integration_workspace_free(w);

    return RelSpeed(v1_x, v2_x);
  };

  double DoubleRarefactionLimit() const {
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    double result, error;
    gsl_function Int;

    // Set up the integral
    struct IntParams Params = {WaveL.S, WaveL.A};
    Int.function = &Integral1;
    Int.params = &Params;

    gsl_integration_qag(&Int, WaveL.p, WaveR.p, 0, 1.0e-12, 1000, 6, w, &result,
                        &error);

    gsl_integration_workspace_free(w);

    return tanh(result);
  };

  double RareShockLimit() const {
    double D = FindD(WaveL, WaveR);
    double h3;
    double J2;
    double Vs;
    double Limit;

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

  void DoubleRarefactionStarValues(double v0) {
    double eps = 1.0e-15;
    double p_min = (WaveR.p + eps) * eps;
    double p_max = WaveL.p;
    double v_star, p_star;
    int status;

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

    do {
      iter++;
      status = gsl_root_fsolver_iterate(s);
      p_star = gsl_root_fsolver_root(s);
      p_min = gsl_root_fsolver_x_lower(s);
      p_max = gsl_root_fsolver_x_upper(s);
      status = gsl_root_test_interval(p_min, p_max, 0, 1.e-12);

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

  void RareShockStarValues(double v12_0) {
    double eps = 1.0e-15;
    double p_min = WaveR.p + eps;
    double p_max = WaveL.p;
    double v_star, p_star;
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

    do {
      iter++;
      status = gsl_root_fsolver_iterate(s);
      p_star = gsl_root_fsolver_root(s);
      p_min = gsl_root_fsolver_x_lower(s);
      p_max = gsl_root_fsolver_x_upper(s);
      status = gsl_root_test_interval(p_min, p_max, 0, 1.e-12);

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

  void TwoShockStarValues(double v12_0) {
    double v_star, p_star;
    double eps = 1.0e-14;
    // double p_min = WaveR.p, p_max = 10.0; // 1.0e11;
    double p_min = WaveL.p + eps, p_max = 1.e12;
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
    gsl_root_fsolver_set(s, &F, p_min, p_max);

    do {
      iter++;
      status = gsl_root_fsolver_iterate(s);
      p_star = gsl_root_fsolver_root(s);
      p_min = gsl_root_fsolver_x_lower(s);
      p_max = gsl_root_fsolver_x_upper(s);
      status = gsl_root_test_interval(p_min, p_max, 0, 1.e-12);

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
    double v12_0 = RelSpeed(WaveL.v, WaveR.v);

    if (v12_0 <= VaccuumLimit()) {
      // This is an unneccessary check (hopefully). Remove when optimizing
      cout << "ERROR: Vaccuum case!" << endl;
      exit(0);

    } else if (v12_0 <= DoubleRarefactionLimit()) {
      // This is the Two Rarefaction case
      double p_star, v_star;
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
  };

  void CalculateIntermediateStates() {
    Wave3.FillWave(&WaveL);
    Wave4.FillWave(&WaveR);
  };

  void FanBoundaries() {
    Wave3.Boundaries(&WaveL, &Wave3, -1.0);
    Wave4.Boundaries(&Wave4, &WaveR, 1.0);
  };
};

int main() {
  cout << setprecision(15);
  RiemannFan Problem;

  // Test Case
  double StateL[4] = {1.0, .0, 0.9, 1000.0};
  double StateR[4] = {1.0, 0.0, 0.9, .01};

  // SR Case
  // double StateL[4] = {1.0, .5, 0.0, 1};
  // double StateR[4] = {.125, 0.0, 0.3, .1};

  // 2R Case
  // double StateL[4] = {1.0, 0.0, 0.9, 1};
  // double StateR[4] = {.125, 0.5, 0.0, .1};

  //  2S Case
  // double StateL[4] = {1.0, 0.5, 0.0, 1};
  // double StateR[4] = {.125, 0.0, 0.999, .1};

  Problem.LoadStates(StateL, StateR);
  Problem.FindWaveTypes();
  Problem.CalculateIntermediateStates();
  Problem.Wave3.PrintWave();
  Problem.Wave4.PrintWave();
  Problem.FanBoundaries();

  cout << "RareFaction Head: " << Problem.Wave3.LeftB * .4 + .5 << endl;
  cout << "RareFaction Tail: " << Problem.Wave3.RightB * .4 + .5 << endl;
  cout << "Contact Point: " << Problem.Wave3.ContactSpeed * .4 + .5 << endl;

  cout << "Shock Location: " << Problem.Wave4.LeftB * .4 + .5 << endl;

  return 0;
}
