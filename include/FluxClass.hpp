#ifndef FLUXCLASS_H_
#define FLUXCLASS_H_

// This is maybe a mistake, but I'm going to define a Flux container to
// hold all of flux data objects and methods

class FluxClass {

public:
  int ndims;
  int Xstart;
  int Xend;
  int Ystart;
  int Yend;
  int nqp;
  bool twoD;
  double *Cons;
  double ***TopFlux;
  double ***BottomFlux;
  double ***RightFlux;
  double ***LeftFlux;
  void Fluxinit(int Ndims, int XStart, int XEnd, int YStart, int YEnd, int Nqp,
                double *COns) {
    ndims = Ndims;
    Xstart = XStart;
    Xend = XEnd;
    Ystart = YStart;
    Yend = YEnd;
    nqp = Nqp;
    Cons = COns;
    if (ndims > 1) {
      twoD = true;

    } else {
      twoD = false;
    }

    LeftFlux = new double **[nqp];
    RightFlux = new double **[nqp];

    for (int i = 0; i < nqp; ++i) {
      LeftFlux[i] = new double *[4];
      RightFlux[i] = new double *[4];
      for (int j = 0; j < 4; ++j) {
        LeftFlux[i][j] = new double[Xend * Yend];
        RightFlux[i][j] = new double[Xend * Yend];
      }
    }

    if (twoD) {

      TopFlux = new double **[nqp];
      BottomFlux = new double **[nqp];

      for (int i = 0; i < nqp; ++i) {
        TopFlux[i] = new double *[4];
        BottomFlux[i] = new double *[4];
        for (int j = 0; j < 4; ++j) {
          TopFlux[i][j] = new double[Xend * Yend];
          BottomFlux[i][j] = new double[Xend * Yend];
        }
      }
    }
  }

  void FOG() {
    for (int qp = 0; qp < nqp; ++qp) {
      for (int var = 0; var < 4; ++var) {
        for (int xdir = 0; xdir < Xend; ++xdir) {
          for (int ydir = 0; ydir < Yend; ++ydir) {
            if (xdir > 0) {
              LeftFlux[qp][var][xdir * Yend + ydir] =
                  Cons[Yend * Xend * var + xdir * Yend + ydir];
            }

            if (xdir < Xend - 1) {
              RightFlux[qp][var][xdir * Yend + ydir] =
                  Cons[Yend * Xend * var + xdir * Yend + ydir];
            }

            if (ydir < Yend - 1) {
              BottomFlux[qp][var][xdir * Yend + ydir] =
                  Cons[Yend * Xend * var + xdir * Yend + ydir];
            }

            if (ydir > 0) {
              TopFlux[qp][var][xdir * Yend + ydir] =
                  Cons[Yend * Xend * var + xdir * Yend + ydir];
            }
          }
        }
      }
    }
  }
};

#endif // FLUXCLASS_H_
