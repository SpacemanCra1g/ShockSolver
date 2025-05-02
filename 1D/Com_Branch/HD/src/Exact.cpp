#include "../include/DomainClass.hpp"
// #include "../include/ExactHLL.hpp"
#include "../include/HydroExact.hpp"
#include <iomanip>

void Domain::Exact(int Start, int Stop) {
  double StateL[3], StateR[3], Result[3];
  double d,p,vx,E;
  
  // std::cout << std::setprecision(15);
  // Prims2Cons(FluxWalls_Prims[LEFT], FluxWalls_Cons[LEFT], Start, Stop + 1);
  // Prims2Cons(FluxWalls_Prims[RIGHT], FluxWalls_Cons[RIGHT], Start, Stop);

  // Prims2Cons(FluxWalls_Prims[LEFT], FluxWalls_Cons[LEFT], Start, Stop);
  // Prims2Cons(FluxWalls_Prims[RIGHT], FluxWalls_Cons[RIGHT], Start, Stop);
  for (int i = Start; i < Stop; ++i) {
    StateL[0] = FluxWalls_Prims[RIGHT][Tidx(DENS, i)];
    StateR[0] = FluxWalls_Prims[LEFT][Tidx(DENS, i + 1)];

    StateL[1] = FluxWalls_Prims[RIGHT][Tidx(VELX, i)];
    StateR[1] = FluxWalls_Prims[LEFT][Tidx(VELX, i + 1)];

    StateL[2] = FluxWalls_Prims[RIGHT][Tidx(PRES, i)];
    StateR[2] = FluxWalls_Prims[LEFT][Tidx(PRES, i + 1)];

    // SolveRiemannFlux(StateL, StateR, Result);
    // if (std::fabs(StateL[3] - StateR[3]) < 1.e-9) {
    //   for (int var = 0; var < 4; ++var) {
    //     Hllc(i, i + 1);
    //   }
    // } else {
      // std::cout << "Cell is  = " << i << std::endl;
      // // if (i == 43) {
      // std::cout << StateL[0] << " " << StateL[1] << " " << StateL[2] << " "
      //           << StateL[3] << " " << std::endl;
      // std::cout << StateR[0] << " " << StateR[1] << " " << StateR[2] << " "
      //           << StateR[3] << " " << std::endl;
      // // }
      // ExactHLL(StateL, StateR, Result);
    ExactSample(StateL, StateR, Result,0.0);

    d = Result[0];
    vx = Result[1];
    p = Result[2];

      // std::cout << "made it passed cell number:  " << i << std::endl;

    CellFlux[Tidx(DENS, i)] = d*vx;
    CellFlux[Tidx(VELX, i)] = d*vx*vx+p;
    E = (d * (0.5 * vx * vx + p / ((GAMMA - 1.0) * d)));
    CellFlux[Tidx(PRES, i)]  = vx*(E + p);

    }
    // Here
  }
