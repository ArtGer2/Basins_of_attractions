#ifndef SYSTEMS_CUH
#define SYSTEMS_CUH


#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <math.h>


#define USE_SYSTEM_FOR_BASINS_2

#ifdef USE_CHAMELEON_MODEL
__device__ inline void calcDiscreteModel(double* X, const double* a, double h) {
    // double h1 = a[0] * h;
    // double h2 = (1 - a[0]) * h;
    // X[0] = X[0] + h1 * (-a[6] * X[1]);
    // X[1] = X[1] + h1 * (a[6] * X[0] + a[1] * X[2]);
    // X[2] = X[2] + h1 * (a[2] - a[3] * X[2] + a[4] * cos(a[5] * X[1]));

    // X[2] = (X[2] + h2 * (a[2] + a[4] * cos(a[5] * X[1]))) / (1 + a[3] * h2);
    // X[1] = X[1] + h2 * (a[6] * X[0] + a[1] * X[2]);
    // X[0] = X[0] + h2 * (-a[6] * X[1]);

    double h1 = a[0] * h;
    double h2 = (1 - a[0]) * h;
    /* Первый этап расчета */
    X[0] = __fma_rn(h1, -a[6] * X[1], X[0]);
    X[1] = __fma_rn(h1, a[6] * X[0] + a[1] * X[2], X[1]);
    double cos_term = cos(a[5] * X[1]);
    X[2] = __fma_rn(h1, a[2] - a[3] * X[2] + a[4] * cos_term, X[2]);

    /* Второй этап расчета */
    X[2] = __fma_rn(h2, (a[2] + a[4] * cos_term), X[2]) / (1 + a[3] * h2);
    X[1] = __fma_rn(h2, (a[6] * X[0] + a[1] * X[2]), X[1]);
    X[0] = __fma_rn(h2, -a[6] * X[1], X[0]);
}
#define SIZE_X 3
#define SIZE_A 7
#define CALC_DISCRETE_MODEL(X, a, h) calcDiscreteModel(X, a, h)
#endif

#ifdef USE_ROSSLER_MODEL
__device__ inline void calcDiscreteModel(double* x, const double* a, double h) {
    double h1 = 0.5 * h + a[0];
    double h2 = 0.5 * h - a[0];

    x[0] = h1 * (-x[1] - x[2]) + x[0];
    x[1] = h1 * (x[0] + a[1] * x[1]) + x[1];
    x[2] = h1 * (a[2] + x[2] * (x[0] - a[3])) + x[2];

    double temp = -h2 * (x[0] - a[3]) + 1.0;
    x[2] = (h2 * a[2] + x[2]) / temp;

    temp = -h2 * a[1] + 1.0;
    x[1] = (h2 * x[0] + x[1]) / temp;

    x[0] = h2 * (-x[1] - x[2]) + x[0];
}
#define SIZE_X 3
#define SIZE_A 4
#define CALC_DISCRETE_MODEL(X, a, h) calcDiscreteModel(X, a, h)
#define CALC_DISCRETE_MODEL_FF(X, a, h) calcDiscreteModelFloatFloat(X, a, h)
#endif



#ifdef USE_SYSTEM_FOR_BASINS
__device__ inline void calcDiscreteModel(double* X, const double* a, double h) {
    float h1 = h * a[0];
    float h2 = h * (1 - a[0]);

    X[0] = X[0] + h * (sin(X[1]) - a[1] * X[0]);
    X[1] = X[1] + h * (sin(X[2]) - a[1] * X[1]);
    X[2] = X[2] + h * (sin(X[0]) - a[1] * X[2]);

    X[2] = (X[2] + h2 * sin(X[0])) / (1 + h2 * a[1]);
    X[1] = (X[1] + h2 * sin(X[2])) / (1 + h2 * a[1]);
    X[0] = (X[0] + h2 * sin(X[1])) / (1 + h2 * a[1]);
}
#define SIZE_X 3
#define SIZE_A 5
#define CALC_DISCRETE_MODEL(X, a, h) calcDiscreteModel(X, a, h)
#endif

#ifdef USE_SYSTEM_FOR_BASINS_2
__device__ inline void calcDiscreteModel(double* X, const double* a, double h) {
    float h1 = h * a[0];
    float h2 = h * (1 - a[0]);

    X[0] = X[0] + h1 * (-X[1]);
    X[1] = X[1] + h1 * (a[1] * X[0] + sin(X[1]));

    float z = X[1];

    X[1] = z + h2 * (a[1] * X[0] + sin(X[1]));
    X[0] = X[0] + h2 * (-X[1]);
}
#define SIZE_X 2
#define SIZE_A 2
#define CALC_DISCRETE_MODEL(X, a, h) calcDiscreteModel(X, a, h)
#endif

#ifdef USE_DAMIR_SYSTEM
__device__ inline void calcDiscreteModel(double* X, const double* a, double h) {
    double X1[3];
    double h1 = __dmul_rn(a[0], h);
    double I;
    double sig;

    sig = 1 - 2 * fmod(floor(2.0 * (X[2]) * a[9]), 2.0);
    I = a[11] + (1.0 * a[10] * (fmod(sig * (2.0 * X[2]), __drcp_rn(a[9]))) * a[9] - 0.5 * sig * a[10]) + 0.5 * a[10];
    //I = a[11] + a[10] * sin(2*3.14159265359 * X[2]*a[9]);
    if (X[0] > a[8]) {
        X[0] = a[3];
        X[1] = __dadd_rn(X[1], a[4]);
    }

    X1[0] = __fma_rn((__fma_rn(__dmul_rn(a[5], X[0]), X[0], __fma_rn(a[6], X[0], a[7])) - X[1] + I), h1, X[0]);
    X1[1] = __fma_rn(__dmul_rn(a[1], __fma_rn(a[2], X[0], -X[1])), h1, X[1]);
    X1[2] = __dadd_rn(X[2], h1);
    //sig = 1 - 2*mod(floor((X1[2])*a[9]),2);
    //I = a[11] + (1.0 * a[12] * (mod(sig*(X1[2]), 1/a[9]))*a[9] - 0.5*sig*a[10]) + 0.5*a[12];
    //I = a[11] + a[10] * sin(2*3.14159265359 * X1[2]*a[9]);
    X[0] = __fma_rn((__fma_rn(__dmul_rn(a[5], X1[0]), X1[0], __fma_rn(a[6], X1[0], a[7])) - X1[1] + I), h, X[0]);
    X[1] = __fma_rn(__dmul_rn(a[1], __fma_rn(a[2], X1[0], -X1[1])), h, X[1]);
    X[2] = __dadd_rn(X[2], h);
}
#define SIZE_X 12
#define SIZE_A 3
#define CALC_DISCRETE_MODEL(X, a, h) calcDiscreteModel(X, a, h)
#endif
#endif // SYSTEMS_CUH
