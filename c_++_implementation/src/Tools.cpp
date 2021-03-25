//
// Created by mateus on 10/03/2021.
//

#include <cmath>
#include <iostream>

#include "../include/Tools.h"

using namespace std;

double
Tools::calculateError(double (*funcAnalytic)(double, double), const double *numericSolution, int nx, const double *domX,
                      double dx,
                      double time, ErrorType errorType) {

    auto *x = new double[nx], *analyticSolution = new double[nx], *difAnalyticToNum = new double[nx];
    double sumNum = 0.0, sumDem = 0.0, maxDifAnalyticToNum = 0.0, maxAnalytic = 0.0;

    for (int i = 0; i < nx; i++) {
        x[i] = domX[0] + i * dx;
        analyticSolution[i] = getAnalyticFuncValue(funcAnalytic, x[i], time);
        difAnalyticToNum[i] = abs(analyticSolution[i] - numericSolution[i]);

        if (difAnalyticToNum[i] > maxDifAnalyticToNum) {
            maxDifAnalyticToNum = difAnalyticToNum[i];
        }

        if (analyticSolution[i] > maxAnalytic) {
            maxAnalytic = analyticSolution[i];
        }

        switch (errorType) {

            case ErrorType::norm1:
                sumNum += abs(analyticSolution[i] - numericSolution[i]);
                sumDem += abs(analyticSolution[i]);
                break;

            case ErrorType::norm2:
                sumNum += pow(analyticSolution[i] - numericSolution[i], 2);
                sumDem += pow(analyticSolution[i], 2);
                break;

            case ErrorType::normInf:
                break;
        }

    }

    switch (errorType) {

        case ErrorType::norm1:
            return sumNum / sumDem;

        case ErrorType::norm2:
            return sqrt(sumNum / sumDem);

        case ErrorType::normInf:
            return maxDifAnalyticToNum / maxAnalytic;

        default:
            return 0.0;
    }
}

double **Tools::readInputFile(const char *fileName) {
    ifstream inputFile;
    inputFile.open(fileName);

    auto **values = new double *[4];

    string aux;
    while (aux != "nPoints") {
        inputFile >> aux;
    }
    int nPoints, readCount = 0, lineCount = 0;
    inputFile >> nPoints;

    for (int i = 1; i < 4; i++) {
        values[i] = new double[nPoints];
    }
    values[0] = new double[1];
    values[0][0] = (double) nPoints;

    double value;
    while (inputFile >> value) {
        if (readCount == 3) {
            readCount = 0;
            lineCount++;
        }
        values[readCount + 1][lineCount] = value;
        readCount++;
    }

    return values;
}