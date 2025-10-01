#ifndef EXACTSOLVER_H
#define EXACTSOLVER_H

typedef double realkind;

void SolveRiemannFlux(realkind StateL[4], realkind StateR[4],
                      realkind Result[4], realkind Time);

#endif
