#include "../include/DomainClass.hpp"
#include "../include/ExactHLL.hpp"
#include "../include/ExactSolver.hpp"
#include <iomanip>

void Domain::Exact(int Start, int Stop) {
  double StateL[4], StateR[4], Result[4];
  double temp;
  double lor;
  double h;
  double val;
  double rho, vx, vy, p;
  std::cout << std::setprecision(15);

  // Prims2Cons(FluxWalls_Prims[LEFT], FluxWalls_Cons[LEFT], Start, Stop);
  // Prims2Cons(FluxWalls_Prims[RIGHT], FluxWalls_Cons[RIGHT], Start, Stop);
  for (int i = Start; i < Stop; ++i) {
    StateL[0] = FluxWalls_Prims[RIGHT][Tidx(DENS, i)];
    StateR[0] = FluxWalls_Prims[LEFT][Tidx(DENS, i + 1)];

    StateL[1] = FluxWalls_Prims[RIGHT][Tidx(VELX, i)];
    StateR[1] = FluxWalls_Prims[LEFT][Tidx(VELX, i + 1)];

    StateL[2] = FluxWalls_Prims[RIGHT][Tidx(VELY, i)];
    StateR[2] = FluxWalls_Prims[LEFT][Tidx(VELY, i + 1)];

    StateL[3] = FluxWalls_Prims[RIGHT][Tidx(PRES, i)];
    StateR[3] = FluxWalls_Prims[LEFT][Tidx(PRES, i + 1)];

    // SolveRiemannFlux(StateL, StateR, Result);
    if (std::fabs(StateL[3] - StateR[3]) < 1.e-9) {
      for (int var = 0; var < 4; ++var) {
        Hllc(i, i + 1);
      }
    } else {
      // std::cout << "Cell is  = " << i << std::endl;
      // // if (i == 43) {
      // std::cout << StateL[0] << " " << StateL[1] << " " << StateL[2] << " "
      //           << StateL[3] << " " << std::endl;
      // std::cout << StateR[0] << " " << StateR[1] << " " << StateR[2] << " "
      //           << StateR[3] << " " << std::endl;
      // // }
      // ExactHLL(StateL, StateR, Result);
      SolveRiemannFlux(StateL, StateR, Result);

      rho = Result[0];
      vx = Result[1];
      vy = Result[2];
      p = Result[3];

      // std::cout << "made it passed cell number:  " << i << std::endl;
      lor = 1.0 / std::sqrt(1.0 - vx * vx - vy * vy);
      h = 1.0 + (GAMMA / (GAMMA - 1.0)) * p / rho;
      val = rho * lor * lor * h;

      CellFlux[Tidx(DENS, i)] = lor * rho * vx;
      CellFlux[Tidx(VELX, i)] = val * vx * vx + p;
      CellFlux[Tidx(VELY, i)] = val * vx * vy;
      CellFlux[Tidx(VELZ, i)] = 0.0;
      CellFlux[Tidx(PRES, i)] = val * vx;
    }
  }
}
