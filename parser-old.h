#pragma once
#ifndef PARSER_OLD_H
#define PARSER_OLD_H
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
static const char *CLEAR_POINT_MATRIX_CMD = "clear";
static const char *SCREEN_RESIZE_CMD = "resize";
static const char *BOX_CMD = "box";
static const char *SPHERE_CMD = "sphere";
static const char *TORUS_CMD = "torus";
static const char *PRINT_TRANSFORMATION_MATRIX_CMD = "p";
static const char *PRINT_POINT_MATRIX_CMD = "pp";
static const char *PLOT_MODE_CMD = "plot-mode";
static const char *CARTESIAN = "cartesian";
static const char *ABSOLUTE = "absolute";
static const char *DRAW_MODE_CMD = "draw-mode";
static const char *LINE_WORD = "line";
static const char *POLYGON_WORD = "polygon";
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
        l: add a line to the edge matrix
            takes 6 arguments (x1, y1, z1, x2, y2, z2)
        i: set the transformation matrix to the identity matrix
        s: create a scale matrix,
            then multiply the transformation matrix by the scale matrix
            takes 3 arguments (sx, sy, sz)
        t: create a translation matrix,
            then multiply the transformation matrix by the translation matrix
            takes 3 arguments (tx, ty, tz)
        x: create an x-axis rotation matrix,
            then multiply the transformation matrix by the rotation matrix
            takes 1 argument (theta)
        y: create an y-axis rotation matrix,
            then multiply the transformation matrix by the rotation matrix
            takes 1 argument (theta)
        z: create an z-axis rotation matrix,
            then multiply the transformation matrix by the rotation matrix
            takes 1 argument (theta)
        a: apply the current transformation matrix to the edge matrix
        c: add a circle to the edge matrix
            takes 4 arguments (cx, cy, cz, r)
        h: add a hermite curve to the edge matrix
            takes 8 arguments (x0, y0, x1, y1, x2, y2, x3, y3)
        b: add a bezier curve to the edge matrix
            takes 8 arguments (x0, y0, x1, y1, x2, y2, x3, y3)
        prism: adds a rectangular prism to the edge matrix
            takes 6 parameters (x, y, z, width, height, depth)
        sphere: adds a sphere to the edge matrix
            takes 5 parameters (x, y, z, radius, step size)
        torus: adds a torus to the edge matrix
            takes 6 parameters (x, y, z, circle_radius, torus_radius, step size)
        p: print the transformation matrix to stdout
        pp: print the point matrix to stdout
        plot-mode: set the plotting mode
            takes 1 parameter (either "cartesian" or "absolute")
        draw-mode: set the drawing mode
            takes 1 parameter (either "line" or "polygon")
        clear: clears the point matrix of all points
        resize: resizes the image/screen and clears the screen
            takes 2 parameters (xres, yres)
        v: draw the lines of the edge matrix to the screen display the screen
        f: save the screen to a file
            takes 1 argument (file name)
        g: draw the lines of the edge matrix to the screen,
            save the screen to a file
            takes 1 argument (file name)
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

/*======== void draw_axes_if_cartesian_mode() ==========
Inputs:
Returns:
Draws the Cartesian axes if the global plotting mode (global_plot_mode) is set
to Cartesian.
=========================================*/
static void draw_axes_if_cartesian_mode();

/*======== void parser_draw() ==========
Inputs:
Returns:
Draws the point matrix to the screen using the global drawing mode
(global_draw_mode).
======================================*/
static void parser_draw();

/*======== void save() ==========
Inputs:     char *line_buf
            int error_is_fatal
Returns:
Saves the generated image to the file specified in the arguments of line_buf.
Returns -1 if no filename was found.
Returns 0 otherwise.
===============================*/
static int save(char *line_buf, int error_is_fatal);

#endif
