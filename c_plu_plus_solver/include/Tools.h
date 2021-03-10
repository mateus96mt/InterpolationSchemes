//
// Created by mateus on 10/03/2021.
//

#ifndef C_PLU_PLUS_SOLVER_TOOLS_H
#define C_PLU_PLUS_SOLVER_TOOLS_H

#include <fstream>
#include <iostream>

using namespace std;

class Tools {

public:

    enum ErrorType {
        norm1, norm2, normInf
    };

    static double getAnalyticFuncValue(double (*funcAnalytic)(double, double), double x, double t) {
        return funcAnalytic(x, t);
    }

    static double
    calculateError(double (*funcAnalytic)(double, double), const double *numericSolution, int nx, const double *domX, double dx,
                   double time, ErrorType errorType);

    static double **readInputFile(const char *fileName);

};

#endif //C_PLU_PLUS_SOLVER_TOOLS_H
