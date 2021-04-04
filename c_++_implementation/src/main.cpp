#include <iostream>
#include <cstring>
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
        if ((x >= 0.0) && (x < 0.2)) {
            return exp(-log(50.0) * ((x - 0.15) / (0.05)) * ((x - 0.15) / (0.05)));
        }
        if ((x > 0.3) && (x < 0.4)) {
            return 1.0;
        }
        if ((x > 0.5) && (x < 0.55)) {
            return 20.0 * x - 10.0;
        }
        if ((x >= 0.55) && (x < 0.6)) {
            return 12.0 - 20.0 * x;
        }
        if ((x > 0.7) && (x < 0.8)) {
            return sqrt(1.0 - ((x - 0.75) / (0.05)) * ((x - 0.75) / (0.05)));
        }
        if ((x > 0.8) && (x <= 2.0)) {
            return 0.0;
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

double initialCond1(double x) {
    return u0(x, 0.0, 1);
}

double initialCond2(double x) {
    return u0(x, 0.0, 2);
}

double initialCond3(double x) {
    return u0(x, 0.0, 3);
}

double analyticSolution1(double x, double t) {
    double a = 1.0;
    return u0(x - a * t, t, 1);
}

double analyticSolution2(double x, double t) {
    double a = 1.0;
    return u0(x - a * t, t, 2);
}

double analyticSolution3(double x, double t) {
    double a = 1.0;
    return u0(x - a * t, t, 3);
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

void testReadFromFile(char *inputFile) {
    auto *solver = new Solver();
    solver->readFile(inputFile);

    switch (solver->getInitialConditionIdentifier()) {
        case 1:
            solver->setFuncInitialCondition(initialCond1);
            break;
        case 2:
            solver->setFuncInitialCondition(initialCond2);
            break;
        case 3:
            solver->setFuncInitialCondition(initialCond3);
            break;
        default:
            solver->setFuncInitialCondition(initialCond1);
            break;
    }

    switch (solver->getAnalyticSolutionIdentifier()) {
        case 1:
            solver->setFuncAnalytic(analyticSolution1);
            break;
        case 2:
            solver->setFuncAnalytic(analyticSolution2);
            break;
        case 3:
            solver->setFuncAnalytic(analyticSolution3);
            break;
        default:
            solver->setFuncAnalytic(analyticSolution1);
            break;
    }

    double (*scheme)(double, double);
    switch (solver->getSchemeIdentifier()) {
        case 1:
            scheme = UpwindSchemes::TOPUS;
            break;
        case 2:
            scheme = UpwindSchemes::FSFL;
            break;
        case 3:
            scheme = UpwindSchemes::EPUS;
            break;
        case 4:
            scheme = UpwindSchemes::SDPUS_C1;
            break;
        case 5:
            scheme = UpwindSchemes::ADBQUICKEST;
            break;
        default:
            scheme = UpwindSchemes::TOPUS;
            break;
    }

    solver->solve(solver->getCfl(), solver->getV(), solver->getA(), "OUTPUT", scheme, 2.0, true);
}

void testReadFromArguments(char *argv[], int c, bool showLog = false) {

    char *output = argv[1];
    EquationType equationType = EquationType::none;
    if (strcmp(argv[2], "linearAdvection") == 0) {
        equationType = EquationType::linearAdvection;
    }
    if (strcmp(argv[2], "boundaryLayer") == 0) {
        equationType = EquationType::boundaryLayer;
    }
    if (strcmp(argv[2], "burges") == 0) {
        equationType = EquationType::burges;
    }
    if (equationType == EquationType::none) {
        printf("Equation type invalid!\n");
        exit(0);
    }

    int nx = stoi(argv[3]);
    auto *domX = new double[2], *domT = new double[2];
    domX[0] = stod(argv[4]);
    domX[1] = stod(argv[5]);
    domT[0] = stod(argv[6]);
    domT[1] = stod(argv[7]);

    double (*_initialCondition)(double);
    switch (stoi(argv[8])) {
        case 1:
            _initialCondition = initialCond1;
            break;
        case 2:
            _initialCondition = initialCond2;
            break;
        case 3:
            _initialCondition = initialCond3;
            break;
        default:
            _initialCondition = initialCond1;
            break;
    }

    double uL, uR;

    if (strcmp(argv[9], "null") == 0) {
        uL = Solver::nan_;
    } else {
        uL = stod(argv[9]);
    }

    if (strcmp(argv[10], "null") == 0) {
        uR = Solver::nan_;
    } else {
        uR = stod(argv[10]);
    }

    double (*_analyticSolution)(double, double);
    switch (stoi(argv[11])) {
        case 1:
            _analyticSolution = analyticSolution1;
            break;
        case 2:
            _analyticSolution = analyticSolution2;
            break;
        case 3:
            _analyticSolution = analyticSolution3;
            break;
        default:
            _analyticSolution = analyticSolution1;
            break;
    }

    double cfl = stod(argv[12]);
    double a = stod(argv[13]);
    double v = stod(argv[14]);
    double (*scheme)(double, double);

    if (strcmp(argv[15], "TOPUS") == 0) {
        scheme = UpwindSchemes::TOPUS;
    }
    if (strcmp(argv[15], "FSFL") == 0) {
        scheme = UpwindSchemes::FSFL;
    }
    if (strcmp(argv[15], "SDPUS_C1") == 0) {
        scheme = UpwindSchemes::SDPUS_C1;
    }
    if (strcmp(argv[15], "EPUS") == 0) {
        scheme = UpwindSchemes::EPUS;
    }
    double param = stod(argv[16]);

    bool isLogActive = c >= 18 && atoi(argv[17]) == 1;

    if (showLog) {
        printf("output: %s\n", output);
        printf("equationType: %d\n", equationType);
        printf("nx: %d\n", nx);
        printf("domX: [%f, %f]\n", domX[0], domX[1]);
        printf("domT: [%f, %f]\n", domT[0], domT[1]);
        printf("initial condition: %d\n", stoi(argv[8]));
        printf("uL: %f\n", uL);
        printf("uR: %f\n", uR);
        printf("analytic solution: %d\n", stoi(argv[11]));
        printf("cfl: %f\n", cfl);
        printf("a: %f\n", a);
        printf("v: %f\n", v);
        printf("scheme: %s\n", argv[15]);
        printf("param: %f\n", param);

    }

    auto *solver = new Solver(nx, domX, domT, _initialCondition, uL, uR,
                              _analyticSolution, equationType);

    solver->solve(cfl, v, a, output, scheme, param, isLogActive);

    double error = Tools::calculateError(_analyticSolution, solver->getValues()[2], nx, domX, solver->getDx(), domT[1],
                                         Tools::ErrorType::norm2);
    printf("%f\n", error);
}

int main(int c, char *argv[]) {


    char *inputFile, *outputFolder;

    testReadFromArguments(argv, c);

//    switch (c) {
//        case 3:
//            inputFile = argv[1];
//            outputFolder = argv[2];
//
//            printf("inputFile: %s\noutputFolder: %s\n", inputFile, outputFolder);
//
//            testReadFromFile(inputFile);
//            break;
//
//        case 10:
//            testReadFromArguments(argv, c);
//            break;
//
//        default:
//            printf("invalid arguments..\n");
//            break;
//    }


    return 0;
}
