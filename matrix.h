#pragma once
#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"

struct matrix {
    double **m;
    int rows, cols;
    int lastcol; // zero-based indexing; index of last empty column
} matrix;

//Basic matrix manipulation routines

/*============== struct matrix *new_matrix() ==============
Inputs:     int rows
            int cols
Returns:
A new matrix with a double **m of dimensions rows x cols

Once allocated, access the matrix as follows:
m->m[r][c]=something;
=========================================================*/

struct matrix *new_matrix(int rows, int cols);

/*============== void free_matrix() ==============
Inputs:     struct matrix *m
Returns:
1. free individual columns
2. free array holding column pointers
3. free actual matrix
=================================================*/
void free_matrix(struct matrix *m);

/*======== void grow_matrix() ==========
Inputs:     struct matrix *m
            int newcols
Returns:
Reallocates the memory for m->m such that it now has
newcols number of columns
======================================*/
void grow_matrix(struct matrix *m, int newcols);

/*============== void copy_matrix() ==============
Inputs:     struct matrix *a
            struct matrix *b
Returns:
Copy matrix a to matrix b
================================================*/
void copy_matrix(struct matrix *a, struct matrix *b);

/*============== void print_matrix() ==============
Inputs:     struct matrix *m
Returns:
Print the matrix
=================================================*/
void print_matrix(struct matrix *m);

/*============== void ident() ==============
Inputs:     struct matrix *m <-- assumes m is a square matrix
Returns:
Converts m in to an identity matrix
==========================================*/
void ident(struct matrix *m);

/*============== void mat_scalar_mult() ==========
Inputs:     double x
            struct matrix *m
Returns:
Multiply each element of m by x
================================================*/
void mat_scalar_mult(double x, struct matrix *m);

/*============== struct matrix *matrix_mult() ==============
Inputs:     struct matrix *a
            struct matrix *b
Returns:
Multiplies matrix a by matrix b
================================================*/
struct matrix *matrix_mult(struct matrix *a, struct matrix *b);

/*====== struct matrix *make_translate() =======
Inputs:     int x
            int y
            int z
Returns:
The translation matrix created using x, y and z
as the translation offsets
================================================*/
struct matrix *make_translate(double x, double y, double z);

/*======== struct matrix *make_scale() ==========
Inputs:     int x
            int y
            int z
Returns:
The dilation matrix created using x, y and z
as the scale factors
================================================*/
struct matrix *make_scale(double x, double y, double z);

/*======== struct matrix *make_rotX() ===========
Inputs:     double theta
Returns:
The rotation matrix created using theta as the
angle of rotation and X as the axis of rotation
================================================*/
struct matrix *make_rotX(double theta);

/*======== struct matrix *make_rotY() ===========
Inputs:     double theta
Returns:
The rotation matrix created using theta as the
angle of rotation and Y as the axis of rotation
================================================*/
struct matrix *make_rotY(double theta);

/*======== struct matrix *make_rotZ() ===========
Inputs:     double theta
Returns:
The rotation matrix created using theta as the
angle of rotation and Z as the axis of rotation.
================================================*/
struct matrix *make_rotZ(double theta);

/*======== struct matrix *apply_transform() ===========
Inputs:     struct matrix *trans_mat
            struct matrix *points
Returns:
The matrix resulting from applying the transformation in trans_mat to the matrix
points.
================================================*/
struct matrix *apply_transform(struct matrix *trans_mat,
                               struct matrix *points);

/*======== struct matrix *make_hermite_coefficients() ===========
Inputs:
Returns:
The coefficient matrix for use in generating cubic Hermite curve equations.
================================================*/
struct matrix *make_hermite_coefficients();

/*======== struct matrix *make_bezier_coefficients() ===========
Inputs:
Returns:
The coefficient matrix for use in generating cubic Bezier curve equations.
================================================*/
struct matrix *make_bezier_coefficients();
#endif
// vim: ts=4:et:sts:sw=4:sr
