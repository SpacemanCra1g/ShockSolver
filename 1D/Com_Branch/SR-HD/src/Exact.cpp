#include "../include/DomainClass.hpp"
#include "../include/ExactHLL.hpp"
#include "../include/QuadExactSolver.hpp"
// #include "../include/ExactSolver.hpp"
#include <iomanip>

void Domain::Exact(int Start, int Stop) {
  realkind StateL[4], StateR[4], Result[4];
  double temp;
  double lor;
  double h;
  double val;
  double rho, vx, vy, p;
  realkind Time = (realkind) 0.0;
  std::cout << std::setprecision(15);

  // Prims2Cons(FluxWalls_Prims[LEFT], FluxWalls_Cons[LEFT], Start, Stop);
  // Prims2Cons(FluxWalls_Prims[RIGHT], FluxWalls_Cons[RIGHT], Start, Stop);
  for (int i = Start; i < Stop; ++i) {
    StateL[0] = (realkind)FluxWalls_Prims[RIGHT][Tidx(DENS, i)];
    StateR[0] = (realkind)FluxWalls_Prims[LEFT][Tidx(DENS, i + 1)];

    StateL[1] = (realkind)FluxWalls_Prims[RIGHT][Tidx(VELX, i)];
    StateR[1] = (realkind)FluxWalls_Prims[LEFT][Tidx(VELX, i + 1)];

    StateL[2] = (realkind)FluxWalls_Prims[RIGHT][Tidx(VELY, i)];
    StateR[2] = (realkind)FluxWalls_Prims[LEFT][Tidx(VELY, i + 1)];

    StateL[3] = (realkind)FluxWalls_Prims[RIGHT][Tidx(PRES, i)];
    StateR[3] = (realkind)FluxWalls_Prims[LEFT][Tidx(PRES, i + 1)];

    // SolveRiemannFlux(StateL, StateR, Result);
    if (std::fabs(StateL[3] - StateR[3]) < 1.e-10) {
        Hllc(i, i + 1);
    } else {
      // std::cout << "Cell is  = " << i << std::endl;
      // // if (i == 43) {
      // std::cout << StateL[0] << " " << StateL[1] << " " << StateL[2] << " "
      //           << StateL[3] << " " << std::endl;
      // std::cout << StateR[0] << " " << StateR[1] << " " << StateR[2] << " "
      //           << StateR[3] << " " << std::endl;
      // // }
      // ExactHLL(StateL, StateR, Result);
      // std::cout << "Cell Number = " << i << std::endl;
      // if (i == 203){
      //   for (int j =0; j < 4; ++j ){
      //     std::cout << StateL[j] << std::endl;
      //   }
      //   std::cout << std::endl;
      //   for (int j =0; j < 4; ++j ){
      //     std::cout << StateR[j] << std::endl;
      //   }
      // }
      SolveRiemannFlux(StateL, StateR, Result, Time);

      rho = (double) Result[0];
      vx = (double) Result[1];
      vy = (double) Result[2];
      p = (double) Result[3];

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
