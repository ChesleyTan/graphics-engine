#pragma once
#ifndef GMATH_H
#define GMATH_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "matrix.h"

/*============== double *calculate_normal() ==============
Inputs:     double x1
            double y1
            double z1
            double x2
            double y2
            double z2
Returns:
A double array of size 3 that represents the vector normal to both
<x1, y1, z1> and <x2, y2, z2>.
========================================================*/
double *calulate_normal(double x1, double y1, double z1,
                        double x2, double y2, double z2);

/*============== double *cross_prod() ==============
Inputs:     double x1
            double y1
            double z1
            double x2
            double y2
            double z2
Returns:
A double array of size 3 that represents the vector cross product of 
<x1, y1, z1> and <x2, y2, z2>.
==================================================*/
double *cross_prod(double x1, double y1, double z1,
                   double x2, double y2, double z2);

/*============== double *dot_prod() ==============
Inputs:     double *u
            double *v
Returns:
The dot product of vectors u and v.
================================================*/
double dot_prod(double *u, double *v);

/*============== double *normalize() =============
Inputs:     double *vector
Returns:
Normalizes the given vector. This operation is performed in-place.
================================================*/
double *normalize(double* vector);

/*============== double *scalar_mult() ===========
Inputs:     double *vector
            double scalar
Returns:
Scales the given vector by the given scalar. This operation is performed in-place.
================================================*/
double *scalar_mult(double* vector, double scalar);

/*============== double *vect_add() ==============
Inputs:     double *u
            double *v
Returns:
Returns a new vector that equals the sum of vectors u and v.
================================================*/
double *vect_add(double* u, double *v);
#endif
