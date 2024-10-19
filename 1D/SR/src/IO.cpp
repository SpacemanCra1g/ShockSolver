#include "../include/DomainClass.hpp"

void Domain::writeResults() {
  FILE *File1 = fopen("OutputData/Density.dat", "w");
  FILE *File2 = fopen("OutputData/VelocityX.dat", "w");
  FILE *File3 = fopen("OutputData/VelocityY.dat", "w");
  FILE *File4 = fopen("OutputData/VelocityZ.dat", "w");
  FILE *File5 = fopen("OutputData/Pressure.dat", "w");
  if (File1 && File2 && File3 && File4 && File5) {

    for (int i = XStart; i < XEnd; i++) {

      fprintf(File1, "%.9g ", DensP[i]);
      fprintf(File2, "%.9g ", Xvel[i]);
      fprintf(File3, "%.9g ", Yvel[i]);
      fprintf(File4, "%.9g ", Zvel[i]);
      fprintf(File5, "%.9g ", Pres[i]);
    }
    fprintf(File1, "\n");
    fprintf(File2, "\n");
    fprintf(File3, "\n");
    fprintf(File4, "\n");
    fprintf(File5, "\n");

    fclose(File1);
    fclose(File2);
    fclose(File3);
    fclose(File4);
    fclose(File5);

  } else {
    printf("There was an issue with the file printing!");
    exit(0);
  }
}
