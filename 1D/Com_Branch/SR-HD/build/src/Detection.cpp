// #include "../include/FluxClass.hpp"

// bool FluxClass::Detection(bool MoodFinished) {

//   double Minlocal, Maxlocal, p;
//   double UL, UR, UM, Mm;
//   int count = 0;

//   bool DMP;

//   double divV, gradP, pL, pR, threshold1, Ca, Mach;

//   MoodFinished = true;

//   // initialize DMP here

//   threshold1 = 5.0;

//   for (int x = 2; x < REdgeX - 2; ++x) {
//     if (MoodOrd[x] == 1) {
//       continue;
//     }

//     if ((Cons[Tidx(Dens, x)] <= 0.) || (std::isnan(Cons[Tidx(Dens, x)]))) {
//       Troubled[x] = true;
//       continue;
//     }

//     p = GetPres(x);

//     if (p <= 0.0 || std::isnan(p)) {
//       Troubled[x] = true;
//       continue;
//     }

//     // First shock detector
//     divV = (Uin[Tidx(MomX, x + 1)] - Uin[Tidx(MomX, x - 1)]) / dx;
//     divV *= 0.5 / Uin[Tidx(Dens, x)];
//     divV -= 0.5 *
//             (Uin[Tidx(MomX, x)] *
//              (Uin[Tidx(Dens, x + 1)] - Uin[Tidx(Dens, x - 1)]) / dx) /
//             std::pow(Uin[Tidx(Dens, x)], 2);

//     // Define Soundspeed
//     Ca = GetPresUin(x) * GAMMA / Uin[(Tidx(Dens, x))];

//     Mach = std::pow(Uin[Tidx(MomX, x)], 2) / pow(Uin[Tidx(Dens, x)], 2);
//     Mach = std::sqrt(Mach / Ca);

//     // Second Shock detector

//     pL = GetPresUin(x - 1);
//     pR = GetPresUin(x + 1);

//     gradP = 0.5 * (std::fabs(pR - pL) / std::fmin(pL, pR)) / dx;
//     if (divV < -dx * dx && Mach > 0.2 && gradP > threshold1) {
//       // std::cout << "DPM Called" << std::endl;
//       // exit(0);
//       DMP = true;
//     } else {
//       DMP = false;
//     }

//     if (DMP) {
//       UL = Uin[Tidx(Dens, x - 1)];
//       UM = Uin[Tidx(Dens, x)];
//       UR = Uin[Tidx(Dens, x + 1)];

//       Minlocal = std::fmin(std::fmin(UL, UM), UR);
//       Maxlocal = std::fmax(std::fmin(UL, UM), UR);

//       Mm = Maxlocal - Minlocal;

//       // Plateau detection
//       if (Mm > dx * dx * dx) {
//         if (Cons[Tidx(Dens, x)] > Maxlocal || Cons[Tidx(Dens, x)] < Minlocal)
//         {
//           Troubled[x] = true;

//           if (false) {
//           }
//         }
//       }
//     }
//   }

//   for (int x = 2; x < REdgeX - 2; ++x) {

//     if (!Troubled[x]) {
//       continue;
//     }
//     MoodFinished = false;

//     MoodOrd[x] -= 2;
//     if (MoodOrd[x] < 0) {
//       std::cout << "1 called trouble" << std::endl;
//       exit(0);
//     }

//     int LOrd, ROrd;
//     LOrd = std::fmin(MoodOrd[(idx(x - 1))], MoodOrd[x]);
//     ROrd = std::fmin(MoodOrd[(idx(x + 1))], MoodOrd[x]);

//     if (LOrd == 3) {
//       for (int var = 0; var < NumVar; ++var) {
//         GPR1Side(Uin, 0, var, x, Left);
//         GPR1Side(Uin, 0, var, x - 1, Right);
//         HLLSide(x - 1);
//       }
//     } else if (LOrd == 1) {
//       for (int var = 0; var < NumVar; ++var) {
//         FOGSide(Uin, 0, var, x, Left);
//         FOGSide(Uin, 0, var, x - 1, Right);
//         HLLSide(x - 1);
//       }
//     }
//     if (ROrd == 3) {
//       for (int var = 0; var < NumVar; ++var) {
//         GPR1Side(Uin, 0, var, x, Right);
//         GPR1Side(Uin, 0, var, x + 1, Left);
//         HLLSide(x + 1);
//       }
//     } else if (ROrd == 1) {
//       for (int var = 0; var < NumVar; ++var) {
//         FOGSide(Uin, 0, var, x, Right);
//         FOGSide(Uin, 0, var, x + 1, Left);
//         HLLSide(x + 1);
//       }
//     }
//   }

//   for (int x = 2; x < REdgeX - 2; ++x) {

//     if (!Troubled[x]) {
//       continue;
//     }
//     for (int var = 0; var < NumVar; ++var) {
//       ReconSide(0, var, x);
//       ReconSide(0, var, x - 1);
//       ReconSide(0, var, x + 1);
//     }
//     Troubled[x] = false;
//     ++count;
//   }
//   std::cout << count << " Troubled Cells" << std::endl;
//   return MoodFinished;
// };
