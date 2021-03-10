//
// Created by mateus on 10/03/2021.
//

#include "../include/Solver.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

Solver::Solver(int nx, double *domX, double *domT, double (*funcInitialCondition)(double), double uL, double uR,
               double (*funcAnalytic)(double, double), EquationType equationType) {

    this->nx = nx;
    this->domX = domX;
    this->domT = domT;
    this->funcInitialCondition = funcInitialCondition;
    this->uL = uL;
    this->uR = uR;
    this->funcAnalytic = funcAnalytic;
    this->equationType = equationType;

    this->numericSolution = new double *[2];
    this->numericSolution[0] = new double[nx];
    this->numericSolution[1] = new double[nx];

    this->xAxes = new double[nx];
    this->analyticSolution = new double[nx];
}

void
Solver::solve(double _cfl, double _v, double _a, const char *_folderName, double (*_schemeFunction)(double, double),
              double _parameter, bool _isLogActive = false) {

    this->cfl = _cfl;
    this->v = _v;
    this->a = _a;
    this->folderName = _folderName;
    this->schemeFunction = _schemeFunction;
    this->parameter = _parameter;
    this->isLogActive = _isLogActive;

    double dx = (domX[1] - domX[0]) / (nx - 1);

    double dt = this->cfl * dx;

    int nt = ((int) ((domT[1] - domT[0]) / dt)) + 1;

    int p = 0, q = 1;

    if (this->isLogActive) {
        this->paramsLog(dx, nt, dt);
    }

    for (int i = 0; i < nx; i++) {
        this->xAxes[i] = domX[0] + i * dx;

        this->numericSolution[p][i] = this->funcInitialCondition(this->xAxes[i]);
        this->numericSolution[q][i] = this->funcInitialCondition(this->xAxes[i]);
    }

    this->numericSolution[p][0] = uL;
    this->numericSolution[q][0] = uL;
    this->numericSolution[p][this->nx - 1] = uR;
    this->numericSolution[q][this->nx - 1] = uR;

    this->saveSolution(this->numericSolution[p], 0, nt, dt, dx);

    for (int t = 1; t < nt; t++) {

        swap(p, q);

        for (int i = 0; i < nx; i++) {

            double velFaceF = 0.0, velFaceG = 0.0, uF = 0.0, uG = 0.0;

            switch (this->equationType) {

                case EquationType::burges:
                    velFaceF = this->meanValueUf(this->numericSolution[q], i);
                    velFaceG = this->meanValueUg(this->numericSolution[q], i);
                    uF = this->meanValueUf(this->numericSolution[q], i);
                    uG = this->meanValueUg(this->numericSolution[q], i);
                    break;

                case EquationType::boundaryLayer:
                case EquationType::linearAdvection:
                    velFaceF = a;
                    velFaceG = a;
                    break;
            }

            double uFaceF = interpolateUf(this->numericSolution[q], i, velFaceF, dx);
            double uFaceG = interpolateUg(this->numericSolution[q], i, velFaceG, dx);

            switch (this->equationType) {

                case EquationType::burges:
                    this->numericSolution[p][i] = ((dt * v *
                                                    (this->numericSolution[q][i - 1] - 2 * this->numericSolution[q][i] +
                                                     this->numericSolution[q][i + 1])) / pow(dx, 2)) -
                                                  (dt * 0.5 * ((uF * uFaceF) - (uG * uFaceG)) / dx) +
                                                  this->numericSolution[q][i];
                    break;

                case EquationType::boundaryLayer:
                    this->numericSolution[p][i] = ((dt * v *
                                                    (this->numericSolution[q][i - 1] - 2 * this->numericSolution[q][i] +
                                                     this->numericSolution[q][i + 1])) / pow(dx, 2)) -
                                                  (dt * a * (uFaceF - uFaceG) / dx) +
                                                  this->numericSolution[q][i];
                    break;

                case EquationType::linearAdvection:

                    this->numericSolution[p][i] = -(dt * a * (uFaceF - uFaceG) / dx) + this->numericSolution[q][i];

            }
        }
        this->saveSolution(this->numericSolution[p], t, nt, dt, dx);
    }

    this->saveSolution(this->numericSolution[p], nt, nt, dt, dx);
    this->solutionIndex = p;

}

double **Solver::getValues() {
    auto **values = new double *[3];
    values[0] = xAxes;
    values[1] = analyticSolution;
    values[2] = numericSolution[solutionIndex];

    return values;
}

double Solver::interpolateUf(const double *u, int i, double velFaceF,
                             double dx) {

    double phiU = 0.0, phiD = 0.0, phiR = 0.0;

    if (velFaceF > 0) {
        phiD = u[i + 1];
        phiR = u[i - 1];
        phiU = u[i];
    } else {
        phiD = u[i];

        if (i + 2 <= nx - 1) {
            phiR = u[i + 2];
        } else {
            phiR = u[i + 1];
        }

        phiU = u[i + 1];
    }

    if ((phiD - phiR) == 0.0 || (phiD - phiR) <= 1e-10) {
        return phiU;
    } else {
        double phiUNorm = normVar(phiU, phiD, phiR);
        double phiFNorm = this->schemeFunction(phiUNorm, this->parameter);
        double phiF = phiR + ((phiD - phiR) * phiFNorm);

        return phiF;
    }

}

double Solver::interpolateUg(const double *u, int i, double velFaceG,
                             double dx) {

    double phiU = 0.0, phiD = 0.0, phiR = 0.0;

    if (velFaceG > 0) {
        phiD = u[i];

        if (i - 2 >= 0) {
            phiR = u[i - 2];
        } else {
            phiR = u[i - 1];
        }

        phiU = u[i - 1];

    } else {
        phiD = u[i - 1];
        phiR = u[i + 1];
        phiU = u[i];
    }

    if ((phiD - phiR) == 0.0 || (phiD - phiR) <= 1e-10) {
        return phiU;
    } else {
        double phiUNorm = normVar(phiU, phiD, phiR);
        double phiFNorm = this->schemeFunction(phiUNorm, this->parameter);
        double phiF = phiR + ((phiD - phiR) * phiFNorm);

        return phiF;
    }

}

void Solver::paramsLog(double dx, int nt, double dt) {
    printf("\n\n------------ PARAMS LOG ------------\n");
    printf("nx: %d\n", nx);
    printf("dx: %f\n", dx);
    printf("nt: %d\n", nt);
    printf("dt: %f\n", dt);
    printf("domX: [%f, %f]\n", this->domX[0], this->domX[1]);
    printf("domT: [%f, %f]\n", this->domT[0], this->domT[1]);
    printf("cfl(param): %f\n", this->cfl);
    printf("cfl(real): %f\n", dt / dx);
    printf("v: %f\n", this->v);
    printf("param: %f\n", this->parameter);
    printf("------------------------------------");
}

void Solver::saveSolution(double *u, int t, int nt, double dt, double dx) {

    double time = this->domT[0] + t * dt;

    //TODO: code to create folder if not exists


    ofstream outputFile;
    char *fileName = new char[100];
    if (t == nt) {
        sprintf(fileName, "%s/FINAL.data", this->folderName);
    } else {
        sprintf(fileName, "%s/%d.data", this->folderName, t);
    }
    outputFile.open(fileName);

    if (this->isLogActive) {
        printf("\n%s", fileName);
    }

    char *params = new char[200];
    sprintf(params, "time = %f    cfl = %f    nx = %d    nt = %d    dx = %f    dt = %f\n", time, this->cfl, this->nx,
            nt, dx, dt);

    outputFile << params;
    outputFile << "nPoints " << this->nx << "\n";

    for (int i = 0; i < this->nx; i++) {
        analyticSolution[i] = funcAnalytic(this->xAxes[i], time);
        outputFile << this->xAxes[i] << " " << this->analyticSolution[i] << " " << u[i] << "\n";
    }

    outputFile.close();
}