#ifndef OPTIONSCLASS_H_
#define OPTIONSCLASS_H_

#include <cstdio>
#include <cstring>
#include <string>

class opt {
public:
  std::string BC, TestProblem;
  int Ndims, ngc, nx, ny, RK, ell, MO, ngp;
  double NegSlope, NN_Thresh, gamma, CFL, TN, T0, y0, yN, x0, xN;
  bool SlowStart, UseNN, UseDMP;

  void ReadInits(std::string FileName) {

    FILE *pFile;
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    char Parameter[50];
    char val[50];
    pFile = fopen(FileName.c_str(), "r");

    while ((read = getline(&line, &len, pFile)) != -1) {
      sscanf(line, "%s %s", Parameter, val);
      if (!strcmp(Parameter, "NDIMS")) {
        Ndims = std::atoi(val);
      }

      else if (!strcmp(Parameter, "NX")) {
        nx = std::atoi(val);
      }

      else if (!strcmp(Parameter, "NY")) {
        ny = std::atoi(val);
      }

      else if (!strcmp(Parameter, "NGC")) {
        ngc = std::atoi(val);
      }

      else if (!strcmp(Parameter, "X0")) {
        x0 = std::atof(val);
      }

      else if (!strcmp(Parameter, "XN")) {
        xN = std::atof(val);
      }

      else if (!strcmp(Parameter, "Y0")) {
        y0 = std::atof(val);
      }

      else if (!strcmp(Parameter, "YN")) {
        yN = std::atof(val);
      }

      else if (!strcmp(Parameter, "T0")) {
        T0 = std::atof(val);
      }

      else if (!strcmp(Parameter, "TN")) {
        TN = std::atof(val);
      }

      else if (!strcmp(Parameter, "CFL")) {
        CFL = std::atof(val);
      }

      else if (!strcmp(Parameter, "gamma")) {
        gamma = std::atof(val);
      }

      else if (!strcmp(Parameter, "ell")) {
        ell = std::atof(val);
      }

      else if (!strcmp(Parameter, "BC")) {
        BC = val;
      }

      else if (!strcmp(Parameter, "TestProblem")) {
        TestProblem = val;
      }

      else if (!strcmp(Parameter, "RK_Method")) {
        RK = std::atoi(val);
      }

      else if (!strcmp(Parameter, "MoodOrder")) {
        MO = std::atoi(val);
      }

      else if (!strcmp(Parameter, "SlowStart")) {
        SlowStart = (!strcmp(val, "True"));
      }

      else if (!strcmp(Parameter, "Use_NN")) {
        UseNN = (!strcmp(val, "True"));
      }

      else if (!strcmp(Parameter, "UseDMP")) {
        UseDMP = (!strcmp(val, "True"));
      }

      else if (!strcmp(Parameter, "NegativeSlope")) {
        NegSlope = std::atof(val);
      }

      else if (!strcmp(Parameter, "NN_Thresh")) {
        NN_Thresh = std::atof(val);
      }

      else if (!strcmp(Parameter, "Num_QuadraturePoints")) {
        ngp = std::atoi(val);
      }

      else if (!strcmp(&Parameter[0], "#")) {
      }

      else {
        printf("There is an issue with the init file, terminating");
        printf("%s \n", Parameter);
        exit(0);
      }
    }
  }
};

#endif // OPTIONSCLASS_H_
