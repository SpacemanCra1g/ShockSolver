#include "../include/DomainClass.hpp"
#include <fstream>

void Domain::WriteOutSolution() {

  std::ofstream file1;

  file1.open("OutputData/Density.dat");
  for (int j = 0; j < REdgeY; ++j) {
    for (int i = 0; i < REdgeX; ++i) {
      file1 << DENS[i * REdgeY + j] << " ";
    }
    file1 << std::endl;
  }
  file1.close();

  file1.open("OutputData/XVel.dat");
  for (int j = 0; j < REdgeY; ++j) {
    for (int i = 0; i < REdgeX; ++i) {
      file1 << XVEL[i * REdgeY + j] << " ";
    }
    file1 << std::endl;
  }
  file1.close();

  file1.open("OutputData/YVel.dat");
  for (int j = 0; j < REdgeY; ++j) {
    for (int i = 0; i < REdgeX; ++i) {
      file1 << YVEL[i * REdgeY + j] << " ";
    }
    file1 << std::endl;
  }
  file1.close();

  file1.open("OutputData/Pres.dat");
  for (int j = 0; j < REdgeY; ++j) {
    for (int i = 0; i < REdgeX; ++i) {
      file1 << PRES[i * REdgeY + j] << " ";
    }
    file1 << std::endl;
  }
  file1.close();
}
