#ifndef CELLCLASS_H_
#define CELLCLASS_H_

#include "GP_Kernel.hpp"
#include <string>

class Cell {
public:
  double *DENS, *PRES, *XVEL, *YVEL, *MOMX, *MOMY, *ENERGY, *gamma;
  double *Cs;
  Cell *LCell, *RCell, *TCell, *BCell;
  int x, y;
  GP_Kernel *GP_Weight;

  // Defined VarConvert.cpp
  void Cons2Prims();
  void Prims2Cons();
  double GetPres();

  void Assign(std::string Var, double val) {
    if (Var == "DENS") {
      *DENS = val;
    } else if (Var == "PRES") {
      *PRES = val;
    } else if (Var == "XVEL") {
      *XVEL = val;
    } else if (Var == "YVEL") {
      *YVEL = val;
    }
  }
};

#endif // CELLCLASS_H_
