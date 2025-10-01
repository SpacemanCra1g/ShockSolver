#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

#define GAMMA (1.4)
#define SIGMA (GAMMA / (GAMMA - 1.0))

using namespace std;

class State {
public:
  double rho, v, p;
  void init(double input[3]) {
    rho = input[0];
    v = input[1];
    p = input[2];
  }
};

class Wave {
public:
  string WaveType;
  double RL, RR, Vx;
  int dir;
  State *Left, *Right;

  void SampleRarefaction(double S, double Out[3]) {
    double a;
    switch (dir) {
    case -1:
      a = sqrt(GAMMA * Left->p / Left->rho);
      Out[0] = Left->rho *
               pow(2.0 / (GAMMA + 1.0) +
                       ((GAMMA - 1.0) / ((GAMMA + 1.0) * a)) * (Left->v - S),
                   2.0 / (GAMMA - 1.0));

      Out[1] = (2.0 / (GAMMA + 1.0)) * (a + .5 * (GAMMA - 1.0) * Left->v + S);

      Out[2] = Left->p *
               pow(2.0 / (GAMMA + 1.0) +
                       ((GAMMA - 1.0) / ((GAMMA + 1.0) * a)) * (Left->v - S),
                   2.0 / (GAMMA - 1.0));
      break;
    case 1:
      a = sqrt(GAMMA * Right->p / Right->rho);
      Out[0] = Right->rho *
               pow(2.0 / (GAMMA + 1.0) -
                       ((GAMMA - 1.0) / ((GAMMA + 1.0) * a)) * (Right->v - S),
                   2.0 / (GAMMA - 1.0));

      Out[1] = (2.0 / (GAMMA + 1.0)) * (-a + .5 * (GAMMA - 1.0) * Right->v + S);

      Out[2] = Right->p *
               pow(2.0 / (GAMMA + 1.0) -
                       ((GAMMA - 1.0) / ((GAMMA + 1.0) * a)) * (Right->v - S),
                   2.0 / (GAMMA - 1.0));
      break;
    };
  };

  void init(State *LeftSide, State *RightSide, int Dir) {
    Left = LeftSide;
    Right = RightSide;
    dir = Dir;

    switch (dir) {
    case -1:
      if (WaveType == "Shock") {
        double AL = 2.0 / ((GAMMA + 1.0) * Left->rho);
        double BL = ((GAMMA - 1.0) / (GAMMA + 1.0)) * Left->p;
        double QL = sqrt((Right->p + BL) / AL);
        Vx = Left->v - QL / Left->rho;
      } else if (WaveType == "Rarefaction") {
        double a = sqrt(GAMMA * Left->p / Left->rho);
        RL = Left->v - a;
        RR = Right->v -
             a * pow(Right->p / Left->p, (GAMMA - 1.0) / (2. * GAMMA));
      } else {
        cout << "Unknown Wavetype, Exiting 67" << endl;
        exit(2);
      }
      break;
    case 1:
      if (WaveType == "Shock") {
        double AR = 2.0 / ((GAMMA + 1.0) * Right->rho);
        double BR = ((GAMMA - 1.0) / (GAMMA + 1.0)) * Right->p;
        double QR = sqrt((Left->p + BR) / AR);
        Vx = Right->v + QR / Right->rho;
      } else if (WaveType == "Rarefaction") {
        double a = sqrt(GAMMA * Right->p / Right->rho);
        RR = Right->v + a;
        RL =
            Left->v + a * pow(Left->p / Right->p, (GAMMA - 1.0) / (2. * GAMMA));
      } else {
        cout << "Unknown Wavetype, Exiting 81" << endl;
        exit(2);
      }
      break;
    default:
      cout << "Unknown wave direction 85" << endl;
      cout << "Exiting" << endl;
      exit(2);
      break;
    };
  };
};

double f(const State &state, const double x) {
  if (state.p < x) {
    double A = (2.0 / ((GAMMA + 1.0) * state.rho));
    double B = ((GAMMA - 1.0) / (GAMMA + 1.0)) * state.p;
    return (x - state.p) * sqrt(A / (x + B));
  } else {
    double a = sqrt(GAMMA * state.p / state.rho);
    a = ((2.0 * a) / (GAMMA - 1.0));
    a *= pow(x / state.p, (GAMMA - 1.0) / (2.0 * GAMMA)) - 1.0;
    return a;
  }
}

double fprime(const State &state, const double x) {
  if (state.p < x) {
    double A = (2.0 / ((GAMMA + 1.0) * state.rho));
    double B = ((GAMMA - 1.0) / (GAMMA + 1.0)) * state.p;
    return sqrt(A / (B + x)) * (1.0 - (x - state.p) / (2.0 * (B + x)));
  } else {
    double a = sqrt(GAMMA * state.p / state.rho);
    a = 1.0 / (state.rho * a);
    a *= pow(x / state.p, -(GAMMA + 1.0) / (2.0 * GAMMA));
    return a;
  }
}

double FindP_Star(const State &Left, const State &Right) {
  double p_star = .5 * (Left.p + Right.p);
  double p_star_new =
      p_star - (f(Left, p_star) + f(Right, p_star) + (Right.v - Left.v)) /
                   (fprime(Left, p_star) + fprime(Right, p_star));

  p_star_new = (p_star_new > 0.0) ? p_star_new : 1.e-10;

  double Difference = fabs(p_star_new - p_star) / (.5 * (p_star + p_star_new));
  int index = 0;

  while (Difference > 1.e-12 && index < 50) {
    p_star = p_star_new;
    p_star_new =
        p_star - (f(Left, p_star) + f(Right, p_star) + (Right.v - Left.v)) /
                     (fprime(Left, p_star) + fprime(Right, p_star));
    Difference = fabs(p_star_new - p_star) / (.5 * (p_star + p_star_new));
    index++;
  }

  if (index == 50) {
    cout << "Newton's Method failed to converge" << endl;
    exit(2);
  }
  return p_star_new;
}

double FindU_Star(State Left, State Right, double p_star) {
  return .5 * (Left.v + Right.v) + .5 * (f(Right, p_star) - f(Left, p_star));
}

void FindRho_Star(State Left, State Right, double p_star, double &Rho_L,
                  double &Rho_R, Wave &LeftWave, Wave &RightWave) {
  double pval, gam;
  if (p_star > Left.p) {
    pval = p_star / Left.p;
    gam = (GAMMA - 1.0) / (GAMMA + 1.0);
    Rho_L = Left.rho * ((pval + gam) / (gam * pval + 1.0));
    LeftWave.WaveType = "Shock";
  } else {
    pval = p_star / Left.p;
    Rho_L = Left.rho * (pow(pval, 1.0 / GAMMA));
    LeftWave.WaveType = "Rarefaction";
  }
  if (p_star > Right.p) {
    pval = p_star / Right.p;
    gam = (GAMMA - 1.0) / (GAMMA + 1.0);
    Rho_R = Right.rho * ((pval + gam) / (gam * pval + 1.0));
    RightWave.WaveType = "Shock";
  } else {
    pval = p_star / Right.p;
    Rho_R = Right.rho * (pow(pval, 1.0 / GAMMA));
    RightWave.WaveType = "Rarefaction";
  }
  return;
}

class RiemannFan {
public:
  State Left, Right, Wave3, Wave4;
  Wave LeftWave, RightWave;

  RiemannFan(double L[3], double R[3]) {
    Left.init(L);
    Right.init(R);

    Wave3.p = FindP_Star(Left, Right);
    Wave3.v = FindU_Star(Left, Right, Wave3.p);

    FindRho_Star(Left, Right, Wave3.p, Wave3.rho, Wave4.rho, LeftWave,
                 RightWave);

    Wave4.p = Wave3.p;
    Wave4.v = Wave3.v;

    LeftWave.init(&Left, &Wave3, -1);
    RightWave.init(&Wave4, &Right, 1);
  };
  void SampleSolution(double S, double Out[3]) {
    if (S < Wave3.v) {
      if (LeftWave.WaveType == "Shock") {
        if (S > LeftWave.Vx) {
          Out[0] = Wave3.rho;
          Out[1] = Wave3.v;
          Out[2] = Wave3.p;
        } else {
          Out[0] = Left.rho;
          Out[1] = Left.v;
          Out[2] = Left.p;
        }
      } else if (LeftWave.WaveType == "Rarefaction") {
        if (S > LeftWave.RR) {
          Out[0] = Wave3.rho;
          Out[1] = Wave3.v;
          Out[2] = Wave3.p;
        } else if (S > LeftWave.RL) {
          LeftWave.SampleRarefaction(S, Out);
        } else {
          Out[0] = Left.rho;
          Out[1] = Left.v;
          Out[2] = Left.p;
        }
      } else {
        cout << "Unknown Wave\n" << "Exit 215" << endl;
        exit(2);
      }
    } else {
      if (RightWave.WaveType == "Shock") {
        if (S < RightWave.Vx) {
          Out[0] = Wave4.rho;
          Out[1] = Wave4.v;
          Out[2] = Wave4.p;
        } else {
          Out[0] = Right.rho;
          Out[1] = Right.v;
          Out[2] = Right.p;
        }
      } else if (RightWave.WaveType == "Rarefaction") {
        if (S < RightWave.RL) {
          Out[0] = Wave4.rho;
          Out[1] = Wave4.v;
          Out[2] = Wave4.p;
        } else if (S < RightWave.RR) {
          RightWave.SampleRarefaction(S, Out);
        } else {
          Out[0] = Right.rho;
          Out[1] = Right.v;
          Out[2] = Right.p;
        }
      } else {
        cout << "Unknown Wave\n" << "Exit 242" << endl;
        exit(2);
      }
    }
  };
};

double VanDerCorput(int i) {
  double result = 0.0;
  double counter = 0.0;
  do {
    result += ((double)(i % 2 != 0)) * pow(.5, counter + 1.0);
    counter++;
    i >>= 1;
  } while (i > 0);
  return result;
}

int main() {
  // int i;
  // for (i = 1; i < 9000 ; ++i){
  //   cout << "i = " << i << " Value is = " << VanDerCorput(i) << endl;
  // }
  double T = .25;
  double L[3] = {1.0, 0.0, 1.};
  double R[3] = {.125, 0.0, 0.1};
  double Dens[400], XVel[400], Pres[400];
  double Out[3];

  RiemannFan Fan(L, R);

  for (int x = 0; x < 400; ++x) {
    Fan.SampleSolution((-.5 + (x / 399.0)) / T, Out);
    Dens[x] = Out[0];
    XVel[x] = Out[1];
    Pres[x] = Out[2];
  }
  FILE *File1 = fopen("OutputData/Density.dat", "w");
  FILE *File2 = fopen("OutputData/VelocityX.dat", "w");
  FILE *File3 = fopen("OutputData/Pressure.dat", "w");

  if (File1 && File2 && File3) {

    for (int i = 0; i < 400; i++) {
      fprintf(File1, "%.9g ", Dens[i]);
      fprintf(File2, "%.9g ", XVel[i]);
      fprintf(File3, "%.9g ", Pres[i]);
    }

    fprintf(File1, "\n");
    fprintf(File2, "\n");
    fprintf(File3, "\n");

    fclose(File1);
    fclose(File2);
    fclose(File3);
  }
  std::cout << 0.25 * Fan.LeftWave.RL + .5 << std::endl;
  std::cout << 0.25 * Fan.LeftWave.RR + .5 << std::endl;
  std::cout << 0.25 * Fan.Wave3.v + .5 << std::endl;
  std::cout << 0.25 * Fan.RightWave.Vx + .5 << std::endl;
  // std::cout << Fan.RightWave.Right << std::endl;

  cout << "Program finished successfully" << endl;
  return 0;
}
