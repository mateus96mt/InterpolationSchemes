#include <iostream>
#include "../include/UpwindSchemes.h"
#include "../include/Tools.h"
#include "../include/Solver.h"

double teste1(double x, double y) {
    return x * y;
}

double teste2(double x, double y) {
    return x + y;
}

double func(double (*f)(double, double), double x, double y) {
    return (*f)(x, y);
}

double funcAnalytic(double x, double t) {
    return pow(x, 2) * t;
}

void testInterpolationSchemes() {
    printf("%f\n", func(teste1, 1, 2));
    printf("%f\n", func(teste2, 1, 2));
    printf("%f\n", func(UpwindSchemes::TOPUS, 0.1, 2.0));
    printf("%f\n", func(UpwindSchemes::FSFL, 0.1, 2.0));
    printf("%f\n", func(UpwindSchemes::SDPUS_C1, 0.1, 2.0));
    printf("%f\n", func(UpwindSchemes::EPUS, 0.1, 2.0));
    printf("%f\n", func(UpwindSchemes::ADBQUICKEST, 0.1, 0.1));
}

void testInputFileRead() {
    double **values = Tools::readInputFile("996.data");

    for (int i = 0; i < (int) values[0][0]; i++) {
        printf("x: %f\ty_ana: %f\ty_num: %f\n", values[1][i], values[2][i], values[3][i]);
    }
}

void testCalculateError() {
    int nx = 100;
    double time = 1.0;
    auto *numeric = new double[nx], *domx = new double[2];
    domx[0] = 0.0;
    domx[1] = 1.0;
    double dx = (domx[1] - domx[0]) / (nx - 1);
    for (int i = 0; i < nx; i++) {
        double x = domx[0] + i * dx;
        numeric[i] = funcAnalytic(x, time) + 0.001;
    }

    printf("Error norm1: %f\n",
           Tools::calculateError(funcAnalytic, numeric, nx, domx, dx, time, Tools::ErrorType::norm1));
    printf("Error norm2: %f\n",
           Tools::calculateError(funcAnalytic, numeric, nx, domx, dx, time, Tools::ErrorType::norm2));
    printf("Error normInf: %f\n",
           Tools::calculateError(funcAnalytic, numeric, nx, domx, dx, time, Tools::ErrorType::normInf));
}

double u0(double x, double t, int _case = 2) {

    if (_case == 1) {
        if (x >= 0 && x < 0.2) {
            double n = pow(((x - 0.15) / 0.05), 2);
            return exp(-(log(n) / log(50)));
        }


        if (x > 0.3 && x < 0.4) {
            return 1.0;
        }

        if (x > 0.5 && x < 0.55) {
            return 20.0 * x - 10.0;
        }

        if (x >= 0.55 && x < 0.66) {
            return -20.0 * x + 12.0;
        }

        if (x > 0.7 && x < 0.8) {
            return sqrt(1 - pow(((x - 0.75) / 0.05), 2));
        }

        return 0.0;
    }

    if (_case == 2) {
        if (x >= 0 && x <= 0.2) {
            return 1.0;
        }

        if (x > 0.2 && x <= 0.4) {
            return 4 * x - 0.6;
        }

        if (x > 0.4 && x <= 0.6) {
            return -4 * x + 2.6;
        }

        if (x > 0.6 && x <= 0.8) {
            return 1.0;
        }

        return 0.0;
    }

    if (_case == 3) {
        if (x >= -1 && x <= -1.0 / 3.0) {
            return -x * sin(3 * M_PI * pow(x, 2) / 2.0);
        }

        if (x > -1.0 / 3.0 && x < 1.0 / 3.0) {
            return abs(sin(2 * M_PI * x));
        }

        if (x >= 1.0 / 3.0 && x <= 1.0) {
            return 2 * x - 1.0 - (1 / 6) * sin(3 * M_PI * x);
        }

        return 0.0;
    }

    return 0.0;
}

double analyticU0(double x, double t, double a, int _case = 2) {
    return u0(x - a * t, _case);
}

double initialCondTest(double x) {
    return u0(x, 0.0, 2);
}

double analyticSolutionTest(double x, double t) {
    double a = 1.0;
    return u0(x - a * t, t, 2);
}

void testSolver() {

    double a = 1.0, v = 0.0, cfl = 0.9;
    auto *domX = new double[2], *domT = new double[2];

    domX[0] = -1.0;
    domX[1] = 1.0;

    domT[0] = 0.0;
    domT[1] = 0.25;

    auto *solver = new Solver(400, domX, domT, initialCondTest, initialCondTest(domX[0]),
                                initialCondTest(domX[1]), analyticSolutionTest, EquationType::linearAdvection);

    solver->solve(cfl, v, a, "OUTPUT", UpwindSchemes::TOPUS, 2.0, true);
}

int main() {

    testSolver();

    return 0;
}
