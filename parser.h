#pragma once
#ifndef PARSER_H
#define PARSER_H
#include <stdio.h>
#include <stdlib.h>

#include "display.h"
#include "draw.h"
#include "matrix.h"

#define LINE_BUF_SIZE 1024
#define CMD_BUF_SIZE 11

// Global variables for use in parse_input() and other parse functions
extern FILE *global_file;
extern struct matrix *global_trans_mat;
extern struct matrix *global_pts;
extern screen global_s;
extern int line_no;

// Use strings for commands for extensibility
static const char *LINE_CMD = "l";
static const char *IDENTITY_CMD = "i";
static const char *SCALE_CMD = "s";
static const char *TRANSLATE_CMD = "t";
static const char *XROT_CMD = "x";
static const char *YROT_CMD = "y";
static const char *ZROT_CMD = "z";
static const char *APPLY_CMD = "a";
static const char *CIRCLE_CMD = "c";
static const char *HERMITE_CURVE_CMD = "h";
static const char *BEZIER_CURVE_CMD = "b";
static const char *VIEW_CMD = "v";
static const char *SAVE_CMD = "f";
static const char *VIEW_AND_SAVE_CMD = "g";
static const char *CLEAR_EDGE_MATRIX_CMD = "clear";
static const char *PRISM_CMD = "prism";
static const char *SPHERE_CMD = "sphere";
static const char *TORUS_CMD = "torus";
static const char *PRINT_TRANSFORMATION_MATRIX_CMD = "p";
static const char *PRINT_POINT_MATRIX_CMD = "pp";
static const char *USE_CARTESIAN_CMD = "cartesian";
static const char *USE_ABSOLUTE_CMD = "absolute";
static const char *QUIT_CMD = "q";

/*======== void parse_input() ==========
Inputs:     char *cmd,
            char *line_buf,
            char error_is_fatal
Returns:
Parses the command according to the following rules:
    Every command occupies its own line
    Any command that requires arguments must have those arguments in the same line.
    The commands are as follows:
        l: add a line to the edge matrix -
            takes 6 arguments (x1, y1, z1, x2, y2, z2)
        i: set the transformation matrix to the identity matrix
        s: create a scale matrix,
            then multiply the transformation matrix by the scale matrix -
            takes 3 arguments (sx, sy, sz)
        t: create a translation matrix,
            then multiply the transformation matrix by the translation matrix -
            takes 3 arguments (tx, ty, tz)
        x: create an x-axis rotation matrix,
            then multiply the transformation matrix by the rotation matrix -
            takes 1 argument (theta)
        y: create an y-axis rotation matrix,
            then multiply the transformation matrix by the rotation matrix -
            takes 1 argument (theta)
        z: create an z-axis rotation matrix,
            then multiply the transformation matrix by the rotation matrix -
            takes 1 argument (theta)
        a: apply the current transformation matrix to the edge matrix
        c: add a circle to the edge matrix -
            takes 4 arguments (cx, cy, cz, r)
        h: add a hermite curve to the edge matrix -
            takes 8 arguments (x0, y0, x1, y1, x2, y2, x3, y3)
        b: add a bezier curve to the edge matrix -
            takes 8 arguments (x0, y0, x1, y1, x2, y2, x3, y3)
        v: draw the lines of the edge matrix to the screen display the screen
        f: save the screen to a file -
            takes 1 argument (file name)
        g: draw the lines of the edge matrix to the screen,
            save the screen to a file -
            takes 1 argument (file name)
        clear: clears the edge matrix of all points
        prism: adds a rectangular prism to the edge matrix -
            takes 6 parameters (x, y, z, width, height, depth)
        sphere: adds a sphere to the edge matrix -
            takes 5 parameters (x, y, z, radius, step size)
        torus: adds a torus to the edge matrix -
            takes 6 parameters (x, y, z, circleRadius, torusRadius, step size)
        p: print the transformation matrix to stdout
        pp: print the point matrix to stdout
        cartesian: set plotting mode to Cartesian
        absolute: set plotting mode to absolute
        q: end parsing (from script file only)

The value of the error_is_fatal argument determines if the program should exit on error.
When running interactively, this may not be the desired behavior.
If the value of error_is_fatal is 0, the program does NOT exit on error.
Otherwise, the program will exit on error.
=====================================*/
void parse_input(char *cmd, char *line_buf, char error_is_fatal);

/*======== void parse_file() ==========
Inputs:   char error_is_fatal
Returns:

Parses file set in synchronize_variables() and performs all of the actions
listed in that file.
The file follows the command format used in parse_input().
The error_is_fatal parameter is passed directly to parse_input().
================================================*/
void parse_file(char error_is_fatal);

/*======== void synchronize_variables() ==========
Inputs:     FILE *file,
            struct matrix *trans_mat,
            struct matrix *pts,
            screen s
Returns:
Hands off the given variables to the parser for use during parsing.
This function is to be called prior to parsing.
=================================================*/
void synchronize_variables(FILE *file,
                           struct matrix *trans_mat,
                           struct matrix *pts,
                           screen s);

/*======== void free_variables() ==========
Inputs:
Returns:
Frees dynamically allocated variables and dereferences pointers used by the parser.
This function is to be called after parsing is completed.
=========================================*/
void free_variables();

#endif
