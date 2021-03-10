#include <iostream>
#include "../include/UpwindSchemes.h"
#include "../include/Tools.h"

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

void testInterpolationSchemes(){
    printf("%f\n", func(teste1, 1, 2));
    printf("%f\n", func(teste2, 1, 2));
    printf("%f\n", func(UpwindSchemes::TOPUS, 0.1, 2.0));
    printf("%f\n", func(UpwindSchemes::FSFL, 0.1, 2.0));
    printf("%f\n", func(UpwindSchemes::SDPUS_C1, 0.1, 2.0));
    printf("%f\n", func(UpwindSchemes::EPUS, 0.1, 2.0));
    printf("%f\n", func(UpwindSchemes::ADBQUICKEST, 0.1, 0.1));
}

void testInputFileRead(){
    double **values = Tools::readInputFile("996.data");

    for (int i = 0; i < (int) values[0][0]; i++) {
        printf("x: %f\ty_ana: %f\ty_num: %f\n", values[1][i], values[2][i], values[3][i]);
    }
}

void testCalculateError(){
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

int main() {






    return 0;
}
