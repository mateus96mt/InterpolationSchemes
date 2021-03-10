//
// Created by mateus on 10/03/2021.
//

#ifndef C_PLU_PLUS_SOLVER_UPWINDSCHEMES_H
#define C_PLU_PLUS_SOLVER_UPWINDSCHEMES_H

#include <math.h>

class UpwindSchemes {

public:

    static double TOPUS(double phiUNorm, double alpha) {
        if (phiUNorm >= 0.0 && phiUNorm <= 1.0) {
            return (alpha * pow(phiUNorm, 4))
                   + (((-2 * alpha) + 1) * pow(phiUNorm, 3))
                   + ((((5 * alpha) - 10) / 4) * pow(phiUNorm, 2))
                   + (((-alpha + 10) / 4) * phiUNorm);
        } else {
            return phiUNorm;
        }
    }

    static double FSFL(double phiUNorm, double beta) {
        if (phiUNorm >= 0.0 && phiUNorm <= 1.0) {
            return (((-2 * beta) + 4) * pow(phiUNorm, 4))
                   + (((4 * beta) - 8) * pow(phiUNorm, 3))
                   + ((((-5 * beta) + 8) / 2) * pow(phiUNorm, 2))
                   + (((beta + 2) / 2) * phiUNorm);
        } else {
            return phiUNorm;
        }
    }

    static double EPUS(double phiUNorm, double lam) {
        if (phiUNorm >= 0.0 && phiUNorm <= 1.0) {
            return (-4 * (lam - 24.0) * pow(phiUNorm, 8))
                   + (16.0 * (lam - 23.0) * pow(phiUNorm, 7))
                   + ((528.0 - (25 * lam)) * pow(phiUNorm, 6))
                   + (((19.0 * lam) - 336.0) * pow(phiUNorm, 5))
                   + ((80.0 - (7.0 * lam)) * pow(phiUNorm, 4))
                   + (lam * pow(phiUNorm, 3))
                   + phiUNorm;
        } else {
            return phiUNorm;
        }
    }

    static double SDPUS_C1(double phiUNorm, double gamma) {
        if (phiUNorm >= 0.0 && phiUNorm <= 1.0) {
            return ((-24.0 + (4 * gamma)) * pow(phiUNorm, 6))
                   + ((68.0 - (12.0 * gamma)) * pow(phiUNorm, 5))
                   + ((-64.0 + (13 * gamma)) * pow(phiUNorm, 4))
                   + ((20.0 - (6 * gamma)) * pow(phiUNorm, 3))
                   + (gamma * pow(phiUNorm, 2))
                   + phiUNorm;
        } else {
            return phiUNorm;
        }
    }

    static double ADBQUICKEST(double phiUNorm, double cfl) {
        double a = (2 - (3 * abs(cfl)) + pow(cfl, 2)) /
                   (7 - (6 * cfl) - (3 * abs(cfl)) + (2 * pow(cfl, 2)));

        double b = (-4 + (6 * cfl) - (3 * abs(cfl)) + pow(cfl, 2)) /
                   (-5 + (6 * cfl) - (3 * abs(cfl)) + (2 * pow(cfl, 2)));

        double phiFNorm = 0;

        if (phiUNorm >= 0 and phiUNorm < a) {
            phiFNorm = (2 - cfl) * phiUNorm;
        }

        if (phiUNorm >= a and phiUNorm <= b) {
            phiFNorm = phiUNorm
                       + 0.5 * (1 - abs(cfl)) * (1 - phiUNorm)
                       - (1 / 6) * (1 - pow(cfl, 2)) * (1 - (2 * phiUNorm));
        }

        if (phiUNorm > b and phiUNorm <= 1) {
            phiFNorm = 1 - cfl + (cfl * phiUNorm);
        }

        if (phiUNorm < 0 or phiUNorm > 1) {
            phiFNorm = phiUNorm;
        }

        return phiFNorm;
    }
};

#endif //C_PLU_PLUS_SOLVER_UPWINDSCHEMES_H
