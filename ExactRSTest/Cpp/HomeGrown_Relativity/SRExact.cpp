#include <cmath>
// #include <gsl/gsl_errno.h>
// #include <gsl/gsl_integration.h>
// #include <gsl/gsl_roots.h>
#include <iomanip>
#include <iostream>
#include <cfenv>


using namespace std;

const long double Gamma = 5.0 / 3;
const long double Sigma = Gamma / (Gamma - 1.0);

long double RelSpeed(long double a, long double b) { return (a - b) / (1 - a * b); }

struct IntParams {
  long double S;
  long double A;
};

struct RarefactionParams {
  struct IntParams Params;
  int sign;
  long double v;
  long double p;
};

struct ShockParams {
  int sign;
  long double h, p, r, v, lor;
};

struct RR_Pstar_Params {
  long double Lv, LA, LS, Lp;
  long double Rv, RA, RS, Rp;
  long double v0;
};

struct RS_Pstar_Params {
  long double v0;
  long double Lv, LA, LS, Lp;
  long double Rh, Rp, Rr, Rlor, Rv;
};
struct SS_Pstar_Params {
  long double Lh, Lp, Lr, Lv, Llor;
  long double Rh, Rp, Rr, Rv, Rlor;
  long double v0;
};
struct RareFactionSample_Params {
  long double xi, pR, pL, SL, AL, vL, sign;
};

long double NewtonIteration(long double (*f)(long double value, void* param),long double x0, long double Tol , void* params){
  long double h, dx;
  h = 1.e-8;
  int count = 0;
  if (x0 - h < 0.0) {x0 = 1.e-3;}
  long double funcEval = f(x0,params);
  while (fabs(funcEval) > Tol && count < 50){
    count++;
  
  // dx = (1.0/12.0)*f(x0-2.*h,params) -(2./3.)*f(x0-h,params) + (2./3.)*f(x0+h,params) - (1./12.)*f(x0+2.*h,params);
  
  dx = -.5*f(x0-h,params) + .5*f(x0+h,params);
  dx /= h;
  x0 -= funcEval/dx;
  // cout << "x0 = " <<x0 <<endl;
  if (x0 - h < 0.0) {x0 = 1.e-3;}
  funcEval = f(x0,params);
  }
  // cout << "count was = " << count << endl;
  if (count == 50) {
    cout << "Newton's method failed" << endl;
  }
  return x0;
};

long double SimpsonsRule(long double (*f)(long double value, void* param),long double min, long double max, long double Tol , void* params){
  long double OldVal;
  long double NewVal;
  long double diff;
  int n = 200000;
  long double h = (max - min)/((long double) n );
  long double runSum = 0.0;
  for (int i = 1; i < n/2; ++i){
    runSum += 4.0 * f(min + (2.0*i - 1.0)*h, params) + 2.0 * f(min + (2.0*i)*h, params);
  }
  runSum += 4.0 * f(min + ((long double) n - 1.0)*h, params);
  runSum += f(min , params) + f(max , params);
  NewVal = runSum*h/3.0;
  return NewVal;

  // do{
  // OldVal = NewVal;
  // n <<= 1;
  // h = (max - min)/((long double) n );
  // runSum = 0.0;
  // for (int i = 1; i < n/2; ++i){
  //   runSum += 4.0 * f(min + (2.0*i - 1.0)*h, params) + 2.0 * f(min + (2.0*i)*h, params);
  // }
  // runSum += 4.0 * f(min + ((long double) n - 1.0)*h, params);
  // runSum += f(min , params) + f(max , params);
  // NewVal = runSum*h/3.0;
  // diff = NewVal - OldVal;
  // } while (fabs(diff) > Tol);
  // cout << "Converged when n = " << n << endl;
  // return OldVal;
};
extern "C" {

static long double Integral1(long double p, void *pram) {
  if (fabs(p) < 1.e-20){
   cout << "P too close to 0 within integral." << endl;
   exit(0);
  }
  struct IntParams *params = (struct IntParams *)pram;
  long double S = params->S;
  long double A = params->A;
  long double rho = pow((p / S), 1.0 / Gamma);
  long double H = 1.0 + Sigma * p / rho;
  long double Cs = sqrt(Gamma * p / (H * rho));
  return sqrt(H * H + A * A * (1 - Cs * Cs)) / ((H * H + A * A) * rho * Cs);
};

long double RarefactionVx(long double p, RarefactionParams *params) {

  // gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
  long double result; //, error;
  // gsl_function Int;

  // Set up the integral
  // Int.function = &Integral1;
  // Int.params = &params->Params;

  // gsl_integration_qag(&Int, params->p, p, 0, 1.0e-12, 1000, 6, w, &result,
  //                     &error);

  // gsl_integration_workspace_free(w);

  result = SimpsonsRule(Integral1 ,params->p, p, 1.e-12 , &params->Params);
  // cout << "131" << endl;

  long double B1 = .5 * log((1.0 + params->v) / (1.0 - params->v));
  return tanh(B1 + (long double)params->sign * result);
};

long double ux(long double xi, long double S, long double press, long double A, long double sign) {
  long double rho = pow(press / S, 1.0 / Gamma);
  long double h = 1.0 + Sigma * press / rho;
  long double cs = sqrt(Gamma * press / (h * rho));
  long double a = cs * h;
  long double b = sign * sqrt(A * A * (1.0 - cs * cs) + h * h);
  return (a - b * xi) / (a * xi - b);
};

long double SampleRarefactionWave(long double pressure, void *Parms) {
  struct RareFactionSample_Params *params =
      (struct RareFactionSample_Params *)Parms;
  struct IntParams IntPar = {params->SL, params->AL};
  struct RarefactionParams RarePar = {IntPar, (int)params->sign, params->vL,
                                      params->pL};

  return ux(params->xi, params->SL, pressure, params->AL, params->sign) -
         RarefactionVx(pressure, &RarePar);
};

long double Taub(long double hA, long double rho, long double p, long double pres) {
  long double c_2 = (1.0 + (p - pres) / (pres * Sigma));
  long double c_1 = -(p - pres) / (pres * Sigma);
  long double c_0 = hA * (p - pres) / rho - hA * hA;
  long double val, error = 1.0;

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

long double J_sqr(long double pres1, long double pres2, long double hA, long double hB) {
  long double val;
  if (fabs(hA - hB) < 1.0e-10) {
    val = Sigma * pres1 * pres2 / (hA * (hA - 1.0));
  } else {
    val = -Sigma * (pres1 - pres2) /
          (hA * (hA - 1.) / pres1 - hB * (hB - 1.0) / pres2);
  }

  if (val == val + 1.0) {
    cout << "landed here" << endl;
    cout << "WaveCrash" << endl;
    long double error = 1 / val;
    error = 1.0 / val;
    error = error / 0.0;
  }
  return val;
};

long double ShockSpeed(long double lor, long double r, long double v, long double J, int sign) {
  long double D = r * lor;
  return (D * D * v + (long double)sign * J * sqrt(J * J + D * D * (1.0 - v * v))) /
         (D * D + J * J);
}

long double ShockVx(long double p, ShockParams *params) {
  long double hA = params->h;
  long double hB = Taub(hA, params->r, params->p, p);
  long double J2 = J_sqr(params->p, p, hA, hB);
  long double J = sqrt(fabs(J2));
  long double Vs = ShockSpeed(params->lor, params->r, params->v, J, params->sign);
  long double Ws = 1.0 / sqrt(1.0 - Vs * Vs);

  return (hA * params->lor * params->v +
          (long double)params->sign * Ws * (p - params->p) / J) /
         (hA * params->lor +
          (p - params->p) * ((long double)params->sign * Ws * params->v / J +
                             1.0 / (params->r * params->lor)));
};

long double DoubleRarefactionPstar(long double p, void *param) {
  struct RR_Pstar_Params *params = (struct RR_Pstar_Params *)param;
  struct IntParams IntPar = {params->LS, params->LA};
  struct RarefactionParams Rare_Vx_Param = {IntPar, -1, params->Lv, params->Lp};

  long double ux3 = RarefactionVx(p, &Rare_Vx_Param);

  IntPar.A = params->RA;
  IntPar.S = params->RS;
  Rare_Vx_Param.Params = IntPar;
  Rare_Vx_Param.sign = 1;
  Rare_Vx_Param.v = params->Rv;
  Rare_Vx_Param.p = params->Rp;

  long double ux4 = RarefactionVx(p, &Rare_Vx_Param);

  long double v13 = RelSpeed(params->Lv, ux3);
  long double v24 = RelSpeed(params->Rv, ux4);

  return RelSpeed(v13, v24) - params->v0;
};

long double RareShockPstar(long double p, void *param) {
  struct RS_Pstar_Params *params = (struct RS_Pstar_Params *)param;
  struct IntParams IntPar = {params->LS, params->LA};
  struct RarefactionParams Rare_Vx_Param = {IntPar, -1, params->Lv, params->Lp};

  long double ux3 = RarefactionVx(p, &Rare_Vx_Param);
  struct ShockParams Shock_Pars;
  Shock_Pars.sign = 1;
  Shock_Pars.h = params->Rh;
  Shock_Pars.p = params->Rp;
  Shock_Pars.r = params->Rr;
  Shock_Pars.v = params->Rv;
  Shock_Pars.lor = params->Rlor;

  long double ux4 = ShockVx(p, &Shock_Pars);

  long double v13 = RelSpeed(params->Lv, ux3);
  long double v64 = RelSpeed(params->Rv, ux4);

  return RelSpeed(v13, v64) - params->v0;
};

long double ShockShockPstar(long double p, void *param) {
  struct SS_Pstar_Params *params = (struct SS_Pstar_Params *)param;

  struct ShockParams LShock = {
      -1, params->Lh, params->Lp, params->Lr, params->Lv, params->Llor,
  };
  struct ShockParams RShock = {
      1, params->Rh, params->Rp, params->Rr, params->Rv, params->Rlor,
  };

  long double ux3 = ShockVx(p, &LShock);
  long double ux4 = ShockVx(p, &RShock);

  long double v13 = RelSpeed(params->Lv, ux3);
  long double v64 = RelSpeed(params->Rv, ux4);
  return RelSpeed(v13, v64) - params->v0;
};
}

class Wave {
public:
  long double rho;
  long double v;
  long double vt;
  long double p;
  string WaveType = "Default";
  long double S;
  long double A;
  long double lor;
  long double h;
  long double ContactSpeed;
  long double LeftB, RightB;
  void SetState(long double state[4]) {
    rho = state[0];
    v = state[1];
    vt = state[2];
    p = state[3];

    SetEntropy();
    Setlor();
    SetHFromState();
    SetA();
  };

  void SampleRare(long double xi, Wave *State, long double sign, long double Sample[4]) {
    long double pmin = fmin(State->p, p);
    long double pmax = fmax(State->p, p);
    long double p_star;
    struct RareFactionSample_Params RareP = {
        xi, p, State->p, State->S, State->A, State->v, sign};
    int status;
    const int max_iter = 150;
    int iter = 0;

    // const gsl_root_fsolver_type *T;
    // gsl_root_fsolver *s;
    // gsl_function F;

    // F.function = &SampleRarefactionWave;
    // F.params = &RareP;
    // T = gsl_root_fsolver_brent;
    // s = gsl_root_fsolver_alloc(T);
    // gsl_root_fsolver_set(s, &F, pmin, pmax);
    p_star = NewtonIteration(SampleRarefactionWave, 0.5*(pmin+pmax), 1.e-12 , &RareP);
    // cout << "319" << endl;

    // do {
    //   iter++;
    //   status = gsl_root_fsolver_iterate(s);
    //   p_star = gsl_root_fsolver_root(s);
    //   pmin = gsl_root_fsolver_x_lower(s);
    //   pmax = gsl_root_fsolver_x_upper(s);
    //   status = gsl_root_test_interval(pmin, pmax, 0, 1.e-12);

    // } while (status == GSL_CONTINUE && iter < max_iter);

    // gsl_root_fsolver_free(s);
    // if (iter == max_iter) {
    //   cout << "Failed to converge before max iteration" << endl;
    //   exit(0);
    // }

    Sample[3] = p_star;
    Sample[0] = pow(p_star / State->S, 1.0 / Gamma);
    long double H = 1.0 + Sigma * p_star / Sample[0];
    // long double ux(long double xi, long double S, long double press, long double A, long double sign) {
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
      long double hB = Taub(Neighbor->h, Neighbor->rho, Neighbor->p, p);
      long double A = Neighbor->A;
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

  long double SVel(Wave *state, long double sS, long double sign) {
    long double SP_sqr = state->v * state->v + state->vt * state->vt;
    long double r = pow(state->p / sS, 1.0 / Gamma);
    long double h = 1.0 + Sigma * state->p / r;
    long double cs = sqrt(Gamma * state->p / (h * r));
    long double cs_sqr = cs * cs;
    return (state->v * (1.0 - cs_sqr) +
            sign * cs *
                sqrt((1.0 - SP_sqr) * (1.0 - SP_sqr * cs_sqr -
                                       state->v * state->v * (1.0 - cs_sqr)))) /
           (1.0 - SP_sqr * cs_sqr);
  };

  void Boundaries(Wave *StateL, Wave *StateR, long double sign) {
    ContactSpeed = v;

    if (WaveType == "Rarefaction") {
      long double s = StateL->p / (pow(StateL->rho, Gamma));
      long double s2 = StateR->p / (pow(StateR->rho, Gamma));
      if (fabs(s - s2) > 1.e-10 * s) {
        cout << "Fail line 256";
        exit(0);
      }

      LeftB = SVel(StateL, s, sign);
      RightB = SVel(StateR, s, sign);

    } else if (WaveType == "Shock") {

      long double hA = StateL->h;
      long double pres = StateR->p;
      long double hB = Taub(hA, StateL->rho, StateL->p, pres);
      long double J2 = J_sqr(StateL->p, StateR->p, hA, hB);
      long double J = sqrt(fabs(J2));
      long double Vs = ShockSpeed(StateL->lor, StateL->rho, StateL->v, J, (int)sign);

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

  void LoadStates(long double StateL[4], long double StateR[4]) {
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
  long double FindD(const Wave &Left, const Wave &Right) const {
    long double D;
    D = 4.0 * Gamma * Left.p;
    D *= ((Gamma - 1.0) * Right.p + Left.p) /
         pow((Gamma - 1.0) * (Left.p - Right.p), 2);
    D *= Right.h * (Right.p - Left.p) / Right.rho - Right.h * Right.h;
    return 1.0 - D;
  };

  long double VaccuumLimit() const {
    long double v1_x, v2_x;
    // gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    long double result; //, error;
    // gsl_function Int;

    // Set up the first integral
    struct IntParams Params = {WaveL.S, WaveL.A};
    // Int.function = &Integral1;
    // Int.params = &Params;

    /*
    gsl_integration_qags is known to work for singularities
    consider replacing the below if encounter difficulties
    */

    //  Assumed no singularities 61pt Gauss-Kronrod
    // gsl_integration_qag(&Int, WaveL.p, 0, 0, 1.0e-13, 1000, 6, w, &result,
                        // &error);
    // gsl_integration_qags(&Int, WaveL.p, 0, 0, 1.0e-13, 1000, w, &result,
    //                      &error);
    result = SimpsonsRule(Integral1, WaveL.p, -0.00001, 1.e-10, &Params);
    // cout << "475" << endl;
    v1_x = tanh(result);

    // Set up the second integral
    Params.A = WaveR.A;
    Params.S = WaveR.S;

    // gsl_integration_qag(&Int, 0, WaveR.p, 0, 1.0e-13, 1000, 6, w, &result,
                        // &error);
    // gsl_integration_qags(&Int, 0, WaveR.p, 0, 1.0e-13, 1000, w, &result,
    //                      &error);
    result = SimpsonsRule(Integral1, 0.00001, WaveR.p, 1.e-10, &Params);
    // cout << "487" << endl;
    v2_x = tanh(result);

    // gsl_integration_workspace_free(w);

    return RelSpeed(v1_x, v2_x);
  };

  long double DoubleRarefactionLimit() const {
    // gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    long double result; //, error;
    // gsl_function Int;

    // Set up the integral
    struct IntParams Params = {WaveL.S, WaveL.A};
    // Int.function = &Integral1;
    // Int.params = &Params;

    // gsl_integration_qag(&Int, WaveL.p, WaveR.p, 0, 1.0e-13, 1000, 6, w, &result,
                        // &error);

    // gsl_integration_workspace_free(w);

    result = SimpsonsRule(Integral1, WaveL.p, WaveR.p, 1.e-10, &Params);
    // cout << "510 " << endl;
    return tanh(result);
  };

  long double RareShockLimit() const {
    long double D = FindD(WaveL, WaveR);
    long double h3;
    long double J2;
    long double Vs;
    long double Limit;

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

  void DoubleRarefactionStarValues(long double v0) {
    long double eps = 1.0e-12;
    long double p_min = (WaveR.p + eps) * eps;
    long double p_max = WaveL.p;
    long double v_star, p_star;
    // int status;

    // cout << "V12_0 Values = " << v0 << endl;
    // if (p_min > p_max) {
    //   cout << "Waves facing the wrong way, 2 Rarefaction case" << endl;
    //   cout << "Terminating" << endl;
    //   exit(0);
    // }

    // const int max_iter = 150;
    // int iter = 0;

    // const gsl_root_fsolver_type *T;
    // gsl_root_fsolver *s;
    // gsl_function F;
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

    // F.function = &DoubleRarefactionPstar;
    // F.params = &Parameters;

    // T = gsl_root_fsolver_brent;
    // s = gsl_root_fsolver_alloc(T);
    // gsl_root_fsolver_set(s, &F, p_min, p_max);
    p_star = NewtonIteration(DoubleRarefactionPstar, .5*(p_min+p_max), 1.e-12, &Parameters);

    // do {
    //   iter++;
    //   status = gsl_root_fsolver_iterate(s);
    //   p_star = gsl_root_fsolver_root(s);
    //   p_min = gsl_root_fsolver_x_lower(s);
    //   p_max = gsl_root_fsolver_x_upper(s);
    //   status = gsl_root_test_interval(p_min, p_max, 0, 1.e-12);

    // } while (status == GSL_CONTINUE && iter < max_iter);

    // gsl_root_fsolver_free(s);
    // if (iter == max_iter) {
    //   cout << "Failed to converge before max iteration" << endl;
    //   exit(0);
    // }

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

  void RareShockStarValues(long double v12_0) {
    long double eps = 1.0e-13;
    long double p_min = WaveR.p - eps;
    long double p_max = WaveL.p;
    long double v_star, p_star;
    // int status;
    // if (p_min > p_max) {
    //   cout << "Waves facing the wrong way, 2 Rarefaction case" << endl;
    //   cout << "Terminating" << endl;
    //   exit(0);
    // }
    // const int max_iter = 150;
    // int iter = 0;
    // const gsl_root_fsolver_type *T;
    // gsl_root_fsolver *s;
    // gsl_function F;
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

    // F.function = &RareShockPstar;
    // F.params = &Parameters;

    // T = gsl_root_fsolver_brent;
    // s = gsl_root_fsolver_alloc(T);
    // gsl_root_fsolver_set(s, &F, p_min, p_max);
    p_star = NewtonIteration(RareShockPstar,.5*(p_min+p_max) , 1.e-13, &Parameters);

    // do {
    //   iter++;
    //   status = gsl_root_fsolver_iterate(s);
    //   p_star = gsl_root_fsolver_root(s);
    //   p_min = gsl_root_fsolver_x_lower(s);
    //   p_max = gsl_root_fsolver_x_upper(s);
    //   status = gsl_root_test_interval(p_min, p_max, 0.0, 1.e-9);

    // } while (status == GSL_CONTINUE && iter < max_iter);

    // gsl_root_fsolver_free(s);
    // if (iter == max_iter) {
    //   cout << "Failed to converge before max iteration" << endl;
    //   exit(0);
    // }
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

  void TwoShockStarValues(long double v12_0) {
    long double v_star, p_star;
    // long double eps = 1.0e-14;
    // long double p_min = WaveR.p, p_max = 10.0; // 1.0e11;
    // long double p_min = WaveL.p - eps, p_max = 1.e12;
    // int status;
    // const int max_iter = 150;
    // int iter = 0;
    // const gsl_root_fsolver_type *T;
    // gsl_root_fsolver *s;
    // gsl_function F;
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

    // F.function = &ShockShockPstar;
    // F.params = &Parameters;

    // T = gsl_root_fsolver_brent;
    // T = gsl_root_fsolver_bisection;
    // T = gsl_root_fsolver_falsepos;
    // s = gsl_root_fsolver_alloc(T);

    // cout << "Min Value = " << ShockShockPstar(p_min, &Parameters) << endl;
    // cout << "Max Value = " << ShockShockPstar(p_max, &Parameters) << endl;
    // exit(0);

    // gsl_root_fsolver_set(s, &F, p_min, p_max);

    p_star = NewtonIteration(ShockShockPstar, 2.*WaveL.p, 1.e-12, &Parameters);

    // do {
    //   iter++;
    //   status = gsl_root_fsolver_iterate(s);
    //   p_star = gsl_root_fsolver_root(s);
    //   p_min = gsl_root_fsolver_x_lower(s);
    //   p_max = gsl_root_fsolver_x_upper(s);
    //   status = gsl_root_test_interval(p_min, p_max, 0, 1.e-12);

    // } while (status == GSL_CONTINUE && iter < max_iter);

    // gsl_root_fsolver_free(s);
    // if (iter == max_iter) {
    //   cout << "Failed to converge before max iteration" << endl;
    //   exit(0);
    // }

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
    long double v12_0 = RelSpeed(WaveL.v, WaveR.v);
    // cout << v12_0 << " The intermediate speed" << endl;
    // cout << "The long double Rare Limit was " << DoubleRarefactionLimit() << endl;
    // cout << "The difference is " << DoubleRarefactionLimit() - v12_0 << endl;
    long double Limit = 0.0; //VaccuumLimit();

    if (v12_0 <= Limit && false) {
      // This is an unneccessary check (hopefully). Remove when optimizing
      cout << "ERROR: Vaccuum case!" << endl;
      exit(0);
    } else {
      Limit = DoubleRarefactionLimit();
      if (v12_0 <= Limit || v12_0 - Limit < 1.e-15) {
        // This is the Two Rarefaction case
        long double p_star, v_star;
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

  void SampleState(long double xi, long double State[4]) {
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

void SolveRiemannFlux(long double StateL[4], long double StateR[4], long double Result[4]) {
  RiemannFan Problem;
  Problem.LoadStates(StateL, StateR);
  Problem.FindWaveTypes();
  Problem.CalculateIntermediateStates();
  Problem.FanBoundaries();
  Problem.SampleState(0.0, Result);
  // delete &Problem;
};

long double quad(long double x, void *par){
  return 3*x*x - 4*x + 12 + log(x)- 10.0;
}

void SolveShockTube(long double StateL[4], long double StateR[4], long double Time) {
  RiemannFan Problem;
  long double Result[4];
  long double Dens[400], XVel[400], YVel[400], Pres[400];
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

int main(){
  // feenableexcept(FE_INVALID);
  // cout << setprecision(15);
  // RiemannFan Problem;

  // Test Case
  long double StateL[4] = {1.0, .0, 0.9, 1000.0};
  long double StateR[4] = {1.0, 0.0, 0.9, .01};
  long double State[4];

  // SolveShockTube(StateL, StateR, .4);
  SolveRiemannFlux(StateL,  StateR, State );
  for (int i = 0; i < 4 ; ++i){
    cout << State[i] << endl;
  }


  // SolveShockTube(StateL, StateR, .4);
  // Problem.LoadStates(StateL, StateR);
  // Problem.FindWaveTypes();
  // Problem.CalculateIntermediateStates();
  // Problem.FanBoundaries();
  // Problem.SampleState(0.0, State);
  // cout << State[3] <<endl;
  
  // cout << Problem.Wave4.rho << endl;
  // cout << Problem.Wave4.v << endl;
  // cout << Problem.Wave4.vt << endl;
  // cout << Problem.Wave4.p << endl;


  // void * blank;
  // long double test = SimpsonsRule(quad,0.5, 4.0, 1.e-12 , blank);
  // cout << test << endl;
  // test = NewtonIteration(quad, 0.10, 1.e-12 , blank);
  // cout << test << endl;
  return 0;


}

