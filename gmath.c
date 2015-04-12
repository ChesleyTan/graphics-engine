#include "gmath.h"

double *calulate_normal(double x1, double y1, double z1,
                        double x2, double y2, double z2) {
    return cross_prod(x1, y1, z1, x2, y2, z2);
}

double *cross_prod(double x1, double y1, double z1,
                   double x2, double y2, double z2) {
    double *cross;
    cross = (double *) malloc(3 * sizeof(double));
    if (cross == NULL) {
        print_error("Memory allocation error.");
    }
    cross[0] = y1*z2 - z1*y2;
    cross[1] = z1*x2 - x1*z2;
    cross[2] = x1*y2 - y1*x2;
    return cross;
}

double dot_prod(double x1, double y1, double z1,
                double x2, double y2, double z2) {
    return x1*x2 + y1*y2 + z1*z2;
}

