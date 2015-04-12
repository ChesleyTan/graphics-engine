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
Inputs:     double x1
            double y1
            double z1
            double x2
            double y2
            double z2
Returns:
The dot product of vectors <x1, y1, z1> and <x2, y2, z2>.
================================================*/
double dot_prod(double x1, double y1, double z1,
                double x2, double y2, double z2);
#endif
