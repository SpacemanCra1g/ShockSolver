#include "../include/IO.hpp"
#include "../include/DomainClass.hpp"

void Domain::writeResults() {
  FILE *File1 = fopen("OutputData/Density.dat", "w");
  FILE *File2 = fopen("OutputData/VelocityX.dat", "w");
#if NDIMS > 1
  FILE *File3 = fopen("OutputData/VelocityY.dat", "w");
#endif
  FILE *File4 = fopen("OutputData/Pressure.dat", "w");
  if (File1 && File2 && File4) {
    /* int EndY; */
    /* if (Mesh->ndims == 1) { */
    /*   EndY = Mesh->StartY + 1; */
    /* } else { */
    /*   EndY = Mesh->EndY; */
    /* } */

    for (int j = YStart; j < YEnd; j++) {
      for (int i = XStart; i < XEnd; i++) {
        /* for(int j = 0; j <Mesh->R_EdgeY;j++){ */
        /*     for(int i =0; i <Mesh->R_EdgeX;i++){ */
        fprintf(File1, "%.9g ", DENS[idx(i, j)]);
        fprintf(File2, "%.9g ", XVEL[idx(i, j)]);
#if NDIMS > 1
        fprintf(File3, "%.9g ", YVEL[idx(i, j)]);
#endif
        fprintf(File4, "%.9g ", PRES[idx(i, j)]);
      }
      fprintf(File1, "\n");
      fprintf(File2, "\n");
#if NDIMS > 1
      fprintf(File3, "\n");
#endif
      fprintf(File4, "\n");
    }
    fclose(File1);
    fclose(File2);
#if NDIMS > 1
    fclose(File3);
#endif
    fclose(File4);
  } else {
    printf("There was an issue with the file printing!");
    exit(0);
  }
}
