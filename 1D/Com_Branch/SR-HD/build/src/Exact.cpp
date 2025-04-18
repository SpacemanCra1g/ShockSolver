#include "../include/DomainClass.hpp"
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

    // if (i == 203) {
    //   std::cout << "Left: " << StateL[0] << " " << StateL[1] << " " <<
    //   StateL[2]
    //             << " " << StateL[3] << " " << std::endl;
    //   std::cout << "Right: " << StateR[0] << " " << StateR[1] << " "
    //             << StateR[2] << " " << StateR[3] << " " << std::endl;
    // }

    if (std::fabs(StateL[3] - StateR[3]) < 1.e-9) {
      Result[0] = 0.5 * (StateL[0] + StateR[0]);
      Result[1] = 0.5 * (StateL[1] + StateR[1]);
      Result[2] = 0.5 * (StateL[2] + StateR[2]);
      Result[3] = 0.5 * (StateL[3] + StateR[3]);

      Result[0] = StateL[0];
      Result[1] = StateL[1];
      Result[2] = StateL[2];
      Result[3] = StateL[3];
      Hllc(i, i + 1);
      rho = Result[0];
      vx = Result[1];
      vy = Result[2];
      p = Result[3];

      // std::cout << "made it passed cell number:  " << i << std::endl;
      lor = 1.0 / std::sqrt(1.0 - vx * vx - vy * vy);
      h = 1.0 + (GAMMA / (GAMMA - 1.0)) * p / rho;
      val = rho * lor * lor * h;

      // CellFlux[Tidx(DENS, i)] = lor * rho * vx;
      // CellFlux[Tidx(VELX, i)] = val * vx * vx + p;
      // CellFlux[Tidx(VELY, i)] = val * vx * vy;
      // CellFlux[Tidx(VELZ, i)] = 0.0;
      // CellFlux[Tidx(PRES, i)] = val * vx;
    } else {
      // std::cout << "I am Cell: " << i << " and I am using the Riemann Solver"
      //           << std::endl;
      // // std::cout << "My difference in pressure is: "
      // //           << std::fabs(StateL[3] - StateR[3]) << std::endl;
      // // std::cout << "Cell Number: " << i << std::endl;
      // std::cout << "Left: " << StateL[0] << " " << StateL[1] << " " <<
      // StateL[2]
      //           << " " << StateL[3] << " " << std::endl;
      // std::cout << "Right: " << StateR[0] << " " << StateR[1] << " "
      //           << StateR[2] << " " << StateR[3] << " " << std::endl;
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
      // std::cout << "Cell = " << i << " P Flux = " << Result[3] << std::endl;

      // Dest[Tidx(DENS, destI)] = lor * d * vx;
      // Dest[Tidx(MOMX, destI)] = val * vx * vx + p;
      // Dest[Tidx(MOMY, destI)] = val * vy * vx;
      // Dest[Tidx(MOMZ, destI)] = val * vz * vx;
      // Dest[Tidx(ENER, destI)] = val * vx;
    }
  }
}
