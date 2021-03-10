//
// Created by mateus on 10/03/2021.
//

#ifndef C_PLU_PLUS_SOLVER_SOLVER_H
#define C_PLU_PLUS_SOLVER_SOLVER_H

enum EquationType {
    burges, boundaryLayer, linearAdvection
};

class Solver {

private:

    double **numericSolution, *analyticSolution, *xAxes;

    int nx, solutionIndex;

    double *domX, *domT, cfl, v, uL, uR, parameter, a, stepInterval = 0.05;

    const char *folderName;

    double (*schemeFunction)(double, double);

    double (*funcAnalytic)(double, double);

    double (*funcInitialCondition)(double);

    bool isLogActive = false;

    EquationType equationType;

    double normVar(double phi, double phiD, double phiR) {
        return (phi - phiR) / (phiD - phiR);
    }

    double meanValueUf(const double *u, int i) {
        return (u[i + 1] + u[i]) / 2.0;
    }

    double meanValueUg(const double *u, int i) {
        return (u[i] + u[-i]) / 2.0;
    }



    double
    interpolateUf(const double *u, int i, double velFaceF,
                  double dx);

    double
    interpolateUg(const double *u, int i, double velFaceG,
                  double dx);

    void paramsLog(double dx, int nt, double dt);

    void saveSolution(double *u, int t, int nt, double dt, double dx);

public:


    Solver(int nx, double *domX, double *domT, double (*funcInitialCondition)(double), double uL, double uR,
           double (*funcAnalytic)(double, double), EquationType equationType);

    void solve(double _cfl, double _v, double _a, const char *_folderName, double (*_schemeFunction)(double, double), double _parameter, bool _isLogActive);

    double **getValues();

};

#endif //C_PLU_PLUS_SOLVER_SOLVER_H
