#include "../include/DomainClass.hpp"

void Domain::PiecewiseLinear(int start, int stop) {
  double Scalefactor, RightWave[NumVar], LeftWave[NumVar];
  double *eigs, **RR;

  // Calculate the sound speed
  Find_Cs(Prims, Cs, start, stop);

  // Calculate the characteristic eigenvectors
  Chars.EigenVectors(Prims, Cs, start, stop);

  // Apply TVD Slope limiting
  (*this.*SlopeLimiter)(Prims, CopyBuffer, &Chars, start, stop);

  // Begin main PLM loop
  for (int i = start; i < stop; ++i) {
    // Pointers to appropriate eigenvectors and eigenvalue arrays
    eigs = Chars.lambda[i];
    RR = Chars.R[i];

    // Reset Right and left wave arrays
    for (int var = 0; var < NumVar; ++var) {
      RightWave[var] = 0.0;
      LeftWave[var] = 0.0;
    }

    for (int lam = 0; lam < NumVar; ++lam) {

      // if lambda_k is moving right
      if (eigs[lam] > 0.0) {
        // Rightwave scale factor
        Scalefactor =
            0.5 * CopyBuffer[Tidx(lam, i)] * (1.0 - eigs[lam] * dt / dx);

        // add the appropriate Right eigenvector * scalefactor
        for (int wave = 0; wave < NumVar; ++wave) {
          RightWave[wave] += Scalefactor * RR[lam][wave];
        }

        // if lambda_k is moving left
      } else {

        // Leftwave scale factor
        Scalefactor =
            0.5 * CopyBuffer[Tidx(lam, i)] * (-1.0 - eigs[lam] * dt / dx);

        // add the appropriate Right eigenvector * scalefactor
        for (int wave = 0; wave < NumVar; ++wave) {
          LeftWave[wave] += Scalefactor * RR[lam][wave];
        }
      }
    }
    // Find the left and rightstate at t = n+1/2
    for (int lam = 0; lam < NumVar; ++lam) {
      FluxWalls_Prims[LEFT][Tidx(lam, i)] = Prims[Tidx(lam, i)] + LeftWave[lam];
      FluxWalls_Prims[RIGHT][Tidx(lam, i)] =
          Prims[Tidx(lam, i)] + RightWave[lam];
    }
  }
}

;
