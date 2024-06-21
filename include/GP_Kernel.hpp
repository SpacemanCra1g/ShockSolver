#ifndef GP_KERNEL_H_
#define GP_KERNEL_H_

class GP_Kernel {
public:
  // Radius 1 prediction vectors
  double *R1_g1xR;
  double *R1_g2xR;

  double *R1_g1xL;
  double *R1_g2xL;

  double *R1_g1yR;
  double *R1_g2yR;

  double *R1_g1yL;
  double *R1_g2yL;

  // Radius 2 prediction vectors
  double *R2_g1xR;
  double *R2_g2xR;
  double *R2_g3xR;

  double *R2_g1xL;
  double *R2_g2xL;
  double *R2_g3xL;

  double *R2_g1yR;
  double *R2_g2yR;
  double *R2_g3yR;

  double R2_g1yL;
  double *R2_g2yL;
  double *R2_g3yL;
};

#endif // GP_KERNEL_H_
