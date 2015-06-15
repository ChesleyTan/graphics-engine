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
        exit(EXIT_FAILURE);
    }
    cross[0] = y1*z2 - z1*y2;
    cross[1] = z1*x2 - x1*z2;
    cross[2] = x1*y2 - y1*x2;
    return cross;
}

double dot_prod(double *u, double*v) {
    return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
}

double *normalize(double* vector) {
    double mag = sqrt(vector[0] * vector[0]
                    + vector[1] * vector[1]
                    + vector[2] * vector[2]);
    if (mag != 0) {
        vector[0] /= mag;
        vector[1] /= mag;
        vector[2] /= mag;
    }
    return vector;
}

double *scalar_mult(double *vector, double scalar) {
    if (isfinite(scalar)) {
        vector[0] *= scalar;
        vector[1] *= scalar;
        vector[2] *= scalar;
    }
    else {
        vector[0] = DBL_MAX;
        vector[1] = DBL_MAX;
        vector[2] = DBL_MAX;
    }
    return vector;
}

double *vect_add(double *u, double *v) {
    double *sum;
    sum = (double *) malloc(3 * sizeof(double));
    if (sum == NULL) {
        print_error("Memory allocation error.");
        exit(EXIT_FAILURE);
    }
    sum[0] = u[0] + v[0];
    sum[1] = u[1] + v[1];
    sum[2] = u[2] + v[2];
    return sum;
}

double *vect_subtract(double *u, double *v) {
    double *difference;
    difference = (double *) malloc(3 * sizeof(double));
    if (difference == NULL) {
        print_error("Memory allocation error.");
        exit(EXIT_FAILURE);
    }
    difference[0] = u[0] - v[0];
    difference[1] = u[1] - v[1];
    difference[2] = u[2] - v[2];
    return difference;
}

double *vect_add_in_place(double *u, double *v) {
    u[0] += v[0];
    u[1] += v[1];
    u[2] += v[2];
    return u;
}

double *clone_vect(double *vector) {
    double *clone;
    clone = (double *) malloc(3 * sizeof(double));
    if (clone == NULL) {
        print_error("Memory allocation error.");
        exit(EXIT_FAILURE);
    }
    clone[0] = vector[0];
    clone[1] = vector[1];
    clone[2] = vector[2];
    return clone;
}

