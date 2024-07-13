#include "../include/DomainClass.hpp"

void Domain::ForwardEuler() {

  Calculate_Quad_Points();

  SolveRiemann();

  TimeStep();
}

void Domain::RK3() {

  SaveDomain();

  DomainAdd(1.0 / 3.0, 2.0 / 3.0);
  ForwardEuler();

  (this->*(this->BC))("Cons");

  ForwardEuler();

  (this->*(this->BC))("Cons");

  DomainAdd(.75, .25);

  ForwardEuler();

  (this->*(this->BC))("Cons");

  DomainAdd(1.0 / 3.0, 2.0 / 3.0);

  Check();
}
