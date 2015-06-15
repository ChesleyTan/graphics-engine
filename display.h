#pragma once
#ifndef DISPLAY_H
#define DISPLAY_H
/*====================== display.h ========================
  Contains functions for basic manipulation of a screen
  represented as a 2 dimensional array of colors.

  A color is an ordered triple of doubles, with each value standing
  for red, green and blue respectively.
=========================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include "utils.h"

#define MAX_COLOR 255

static int XRES = 500;
static int YRES = 500;
/*
 * XRES_CARTESIAN and YRES_CARTESIAN should be half of XRES AND YRES
 * respectively.
 */
static int XRES_CARTESIAN = 250;
static int YRES_CARTESIAN = 250;

extern double **z_buffer;

/*
 * Every point has an individual int for each color value
 */
struct point_t {
    double red;
    double green;
    double blue;
} point_t;

typedef struct point_t color;

// Used for initializing new color structs to zero
static const color init_color;

typedef struct point_t **screen;

/*======== screen resize_screen() ==========
Inputs:     screen s,
            int x_res,
            int y_res
Returns:
Sets the XRES and YRES variables to x_res and y_res respectively.
XRES and YRES determines the screen size.
Takes a screen, frees it, and returns a new screen with the new dimensions.
NOTE: All points already drawn in the screen will be lost!
==========================================*/
screen resize_screen(screen s, int x_res, int y_res);

/*======== screen new_screen() ==========
Inputs:
Returns:
A new screen of size XRES * YRES.
=======================================*/
screen new_screen();

/*======== void clear_screen() ==========
Inputs:     screen s
Returns:
Sets every color in screen s to black.
=======================================*/
void clear_screen(screen s);

/*======== void free_screen() ==========
Inputs:     screen s
Returns:
Frees all pointers in screen s.
=======================================*/
void free_screen(screen s);

/*======== void plot_absolute() ==========
Inputs:     screen s
            color c
            int x
            int y
            int z
Returns:
Sets the color at pixel x, y to the color represented by c.
NOTE: s[0][0] will be the lower left hand corner of the screen.
========================================*/
void plot_absolute(screen s, color c, int x, int y, int z);

/*======== void plot_cartesian() ==========
Inputs:     screen s
            color c
            int x
            int y
            int z
Returns:
Sets the color at pixel x, y in the standard cartesian plane
to the color represented by c.
NOTE: s[0][0] will be the center of the screen.
========================================*/
void plot_cartesian(screen s, color c, int x, int y, int z);

/*======== void plot_cartesian_wrap() ==========
Inputs:     screen s
            color c
            int x
            int y
            int z
Returns:
Sets the color at pixel x, y in the standard cartesian plane
to the color represented by c.
NOTE: s[0][0] will be the center of the screen.
Coordinates that exceed the dimensions of the screen array will be automatically
wrapped.
========================================*/
void plot_cartesian_wrap(screen s, color c, int x, int y, int z);

/*======== void save_ppm() ==========
Inputs:     screen s
            char *file
Returns:
Saves screen s as a valid ppm file using the
settings in display.h.
===================================*/
void save_ppm(screen s, char *file);

/*======== void save_extension() ==========
Inputs:     screen s
            char *file
Returns:
Saves the screen stored in s to the filename represented
by file.
If the extension for file is an image format supported
by the "convert" command, the image will be saved in
that format.
=========================================*/
void save_extension(screen s, char *file);

/*======== void display() ==========
Inputs:     screen s
Returns:
Displays the screen s.
==================================*/
void display(screen s);

/*======== void allocate_z_buffer() ==========
Inputs:
Returns:
Allocates memory for the Z buffer.
============================================*/
void allocate_z_buffer();

/*======== void free_z_buffer() ==========
Inputs:
Returns:
Frees memory allocated for the Z buffer.
========================================*/
void free_z_buffer();

/*======== color avg_color() ==========
Inputs:     color c1
            color c2
Returns:
The average color between c1 and c2 (arithmetic mean of individual color
components).
======================================*/
color avg_color(color c1, color c2);

/*======== color add_color() ==========
Inputs:     color c1
            color c2
Returns:
The sum of colors c1 and c2.
=====================================*/
color add_color(color c1, color c2);

/*======== color subtract_colors() ==========
Inputs:     color c1
            color c2
Returns:
The difference of colors c1 and c2.
===========================================*/
color subtract_color(color c1, color c2);

/*======== color divide_color() ==========
Inputs:     color c
            int n
Returns:
Divides each color component of c by n.
========================================*/
color divide_color(color c, int n);

/*======== void constrain_color() ==========
Inputs:     color *c
Returns:
Limits the given color to values between 0 and MAX_COLOR, inclusive.
This is performed in-place.
==========================================*/
void constrain_color(color *c);

/*======== void swap_colors() ==========
Inputs:     color *c1
            color *c2
Returns:
Swaps the corresponding color components for c1 and c2.
==========================================*/
void swap_colors(color *c1, color *c2);
#endif
// vim: ts=4:et:sts:sw=4:sr
