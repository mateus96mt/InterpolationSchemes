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
    this->funcInterpolationScheme = _schemeFunction;
    this->parameter = _parameter;
    this->isLogActive = _isLogActive;

    double dx = (domX[1] - domX[0]) / (nx - 1);
    this->dx = dx;

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

    if (uL != nan_) {
        this->numericSolution[p][0] = uL;
        this->numericSolution[q][0] = uL;
    } else {
        this->numericSolution[p][0] = funcInitialCondition(this->xAxes[0]);
        this->numericSolution[q][0] = funcInitialCondition(this->xAxes[0]);
    }


    if (uR != nan_) {
        this->numericSolution[p][this->nx - 1] = uR;
        this->numericSolution[q][this->nx - 1] = uR;
    } else {
        this->numericSolution[p][this->nx - 1] = funcInitialCondition(this->xAxes[this->nx - 1]);
        this->numericSolution[q][this->nx - 1] = funcInitialCondition(this->xAxes[this->nx - 1]);
    }


    this->saveSolution(this->numericSolution[p], 0, nt, dt, dx);

    for (int t = 1; t < nt; t++) {

        swap(p, q);

        for (int i = 1; i < nx-1; i++) {

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
        if (phiUNorm < 0.0 || phiUNorm > 1.0) {
            return phiU;
        } else {
            double phiFNorm = this->funcInterpolationScheme(phiUNorm, this->parameter);
            double phiF = phiR + ((phiD - phiR) * phiFNorm);

            return phiF;
        }

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
        if (phiUNorm < 0.0 || phiUNorm > 1.0) {
            return phiU;
        } else {
            double phiFNorm = this->funcInterpolationScheme(phiUNorm, this->parameter);
            double phiF = phiR + ((phiD - phiR) * phiFNorm);

            return phiF;
        }
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
    printf("------------------------------------\n\n");
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
        printf("%s\n", fileName);
    }

    char *params = new char[200];
    sprintf(params, "time = %f    cfl = %f    nx = %d    nt = %d    dx = %f    dt = %f\n", time, this->cfl, this->nx,
            nt, dx, dt);

    outputFile << params;
    outputFile << "nPoints " << this->nx << "\n";

    for (int i = 0; i < this->nx; i++) {
        if (funcAnalytic != nullptr) {
            analyticSolution[i] = funcAnalytic(this->xAxes[i], time);
        } else {
            analyticSolution[i] = 0.0;
        }

        outputFile << this->xAxes[i] << " " << this->analyticSolution[i] << " " << u[i] << "\n";
    }

    outputFile.close();
}

void Solver::readFile(char *fileName) {
    try {
        ifstream inputFile;
        inputFile.open(fileName);

        string aux;
        inputFile >> aux >> aux;
        if (aux == "linearAdvection") {
            this->equationType = EquationType::linearAdvection;
        }
        if (aux == "boundaryLayer") {
            this->equationType = EquationType::boundaryLayer;
        }
        if (aux == "burges") {
            this->equationType = EquationType::burges;
        }

        inputFile >> aux >> this->nx;
        this->numericSolution = new double *[2];
        this->numericSolution[0] = new double[this->nx];
        this->numericSolution[1] = new double[this->nx];

        this->xAxes = new double[this->nx];
        this->analyticSolution = new double[this->nx];

        this->domX = new double[2];
        inputFile >> aux >> this->domX[0] >> this->domX[1];

        this->domT = new double[2];
        inputFile >> aux >> this->domT[0] >> this->domT[1];

        inputFile >> aux >> this->initialConditionIdentifier;

        inputFile >> aux >> aux;
        if (aux == "null") {
            this->uL = nan_;
        }

        inputFile >> aux >> aux;
        if (aux == "null") {
            this->uR = nan_;
        }

        inputFile >> aux >> this->analyticSolutionIdentifier;

        inputFile >> aux >> this->cfl;
        inputFile >> aux >> this->a;
        inputFile >> aux >> this->v;

        inputFile >> aux >> aux;
        cout << aux;
        if (aux == "TOPUS") {
            this->schemeIdentifier = 1;
        }
        if (aux == "FSFL") {
            this->schemeIdentifier = 2;
        }
        if (aux == "EPUS") {
            this->schemeIdentifier = 3;
        }
        if (aux == "SDPUS_C1") {
            this->schemeIdentifier = 4;
        }
        if (aux == "ADBQUICKEST") {
            this->schemeIdentifier = 5;
        }

        inputFile >> this->parameter;

        printf("\ninitialConditionIdentifier: %d", this->initialConditionIdentifier);
        printf("\nanalyticSolutionIdentifier: %d", this->analyticSolutionIdentifier);
        printf("\nschemeIdentifier: %d", this->schemeIdentifier);


    } catch (const exception e) {
        printf("\nError in input file %s, example of valid input (may be outdated):\n", fileName);
        printf("\nEquationType linearAdvection\ndomX -1.0 1.0\ndomT 0.0 0.25\ninitialCond 1\nuL null\nuR null\nanalyticSolution 1\ncfl 0.9\na 1.0\nv 0.0\ninterpolationScheme TOPUS 2.0\n");
    }
}

int Solver::getInitialConditionIdentifier() const {
    return initialConditionIdentifier;
}

int Solver::getAnalyticSolutionIdentifier() const {
    return analyticSolutionIdentifier;
}

int Solver::getSchemeIdentifier() const {
    return schemeIdentifier;
}

void Solver::setFuncAnalytic(double (*_funcAnalytic)(double, double)) {
    Solver::funcAnalytic = _funcAnalytic;
}

void Solver::setFuncInitialCondition(double (*_funcInitialCondition)(double)) {
    Solver::funcInitialCondition = _funcInitialCondition;
}

double Solver::getCfl() const {
    return cfl;
}

double Solver::getV() const {
    return v;
}

double Solver::getA() const {
    return a;
}

double Solver::getDx() const {
    return dx;
}
