#pragma once
#ifndef DRAW_H
#define DRAW_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "display.h"
#include "gmath.h"
#include "matrix.h"
#include "utils.h"

/* Step size for drawing circles and curves */
#define STEP_SIZE 0.001
#define MIN_STEP_SIZE 0.0001

// TODO doc
extern double view_vector[3];
extern double light_vector[3];
extern int i_ambient_r;
extern int i_ambient_g;
extern int i_ambient_b;
extern double k_ambient_r;
extern double k_ambient_g;
extern double k_ambient_b;
extern int i_point_r;
extern int i_point_g;
extern int i_point_b;
extern double k_point_r;
extern double k_point_g;
extern double k_point_b;

/* enum for plotting mode for use in functions defined in draw.c
 * Its value can be one of PLOT_CARTESIAN or PLOT_ABSOLUTE.
 * If PLOT_CARTESIAN is selected, the points are plotted on the Cartesian axes,
 * with the origin at the center.
 * If PLOT_ABSOLUTE is selected, the points are plotted with the origin at the
 * bottom-left corner.
 * */
typedef enum {
    PLOT_CARTESIAN,
    PLOT_ABSOLUTE
} plotting_mode;

/* enum for drawing mode for use in functions defined in draw.c
 * Its value can be one of DRAW_LINE or DRAW_POLYGON.
 * If DRAW_LINE is selected, points will be paired to form lines/edges.
 * If DRAW_POLYGON is selected, points will be grouped in triples to form
 * triangles.
 * */
typedef enum {
    DRAW_LINE,
    DRAW_POLYGON
} drawing_mode;

/* enum for rendering mode for use in functions defined in draw.c
 * Its value can be one of RENDER_WIREFRAME or RENDER_SURFACE.
 * If RENDER_WIREFRAME is selected, polygons will be rendered as wireframes.
 * If RENDER_SURFACE is selected, points will be rendered as surfaces.
 */
typedef enum {
    RENDER_WIREFRAME,
    RENDER_SURFACE
} rendering_mode;

/* enum for curve types passed into add_curve() */
typedef enum {
    HERMITE_CURVE,
    BEZIER_CURVE
} curve_type;

/* Plotting mode to be used by default globally.
 * This value may be set programmatically. */
extern plotting_mode global_plot_mode;

/* Drawing mode to be used by default globally.
 * This value may be set programmatically. */
extern drawing_mode global_draw_mode;

/* Rendering mode to be used by default globally.
 * This value may be set programmatically. */
extern rendering_mode global_render_mode;

/*======== void add_point() ==========
Inputs:     struct matrix *points
            double x
            double y
            double z
Returns:
Adds point (x, y, z) to 'points' and increments 'points.lastcol'.
If 'points' is full, 'points' is automatically resized.
====================================*/
void add_point(struct matrix * points, double x, double y, double z);

/*======== void add_edge() ==========
Inputs:     struct matrix *points
            double x0, double y0, double z0, double x1, double y1, double z1
Returns:
Add the line connecting (x0, y0, z0) to (x1, y1, z1) to 'points'.
===================================*/
void add_edge(struct matrix * points,
            double x0, double y0, double z0,
            double x1, double y1, double z1);

/*======== void add_circle() ==========
Inputs:     struct matrix *points
            double cx, double cy, double cz, double r, double step
Returns:
Add the circle centered at (cx, cy, cz) with radius r to 'points', using the step
size defined by the step parameter.
===================================*/
void add_circle(struct matrix *points,
                double cx,
                double cy,
                double cz,
                double r,
                double step);

/*======== void add_curve() ==========
Inputs:     struct matrix *points
            double step,
            curve_type type,
            double x0, double y0,
            double x1, double y1,
            double x2, double y2,
            double x3, double y3
Returns:
Adds a curve of type defined by the type parameter to 'points', using a step size
defined by the step parameter, and control points (x0, y0), (x1, y1), (x2, y2),
(x3, y3). The usage of these points differs based on the curve type chosen. For
a list of curve types, see the curve_type enum.
====================================*/
void add_curve(struct matrix *points,
               double step,
               curve_type type,
               double x0, double y0,
               double x1, double y1,
               double x2, double y2,
               double x3, double y3);

/*======== void add_hermite_curve() ==========
Inputs:     struct matrix *points
            double step,
            double x0, double y0,
            double x1, double y1,
            double dx0, double dy0,
            double dx1, double dy1
Returns:
Adds a cubic Hermite curve to 'points', using a step size defined by the
step parameter, endpoints (x0, y0) and (x1, y1), and first
derivatives at the endpoints (dx0, dy0) and (dx1, dy1).
====================================*/
void add_hermite_curve(struct matrix *points,
                       double step,
                       double x0, double y0,
                       double x1, double y1,
                       double dx0, double dy0,
                       double dx1, double dy1);

/*======== void add_bezier_curve() ==========
Inputs:     struct matrix *points
            double step,
            double x0, double y0,
            double x1, double y1,
            double x2, double y2,
            double x3, double y3
Returns:
Adds a cubic Bezier curve to 'points', using a step size defined by the
step parameter and control points (x0, y0), (x1, y1), (x2, y2), and (x3, y3).
====================================*/
void add_bezier_curve(struct matrix *points,
                      double step,
                      double x0, double y0,
                      double x1, double y1,
                      double x2, double y2,
                      double x3, double y3);

/*======== void add_box() ==========
Inputs:     struct matrix *points,
            double x,
            double y,
            double z,
            double width,
            double height,
            double depth,
            drawing_mode draw_mode
Returns:
Adds the corners of a rectangular prism with upper-left corner of its front face
at (x, y, z) and width, height, and depth equal to those given to the
matrix 'points'. The draw_mode parameter determines the method used for drawing
the prism. See drawing_mode for more information.
==================================*/
void add_box(struct matrix *points,
             double x,
             double y,
             double z,
             double width,
             double height,
             double depth,
             drawing_mode draw_mode);

/*======== void add_sphere() ==========
Inputs:     struct matrix *points,
            double step,
            double x,
            double y,
            double z,
            double radius,
            drawing_mode draw_mode
Returns:
Adds the edges/polygons of a sphere centered at (x, y, z) with a radius equal to
the one given to the matrix 'points' using the given step size 
(0 < step < 2).
The draw_mode parameter determines the method used for drawing the prism. See
drawing_mode for more information.
====================================*/
void add_sphere(struct matrix *points,
                double step,
                double x,
                double y,
                double z,
                double radius,
                drawing_mode draw_mode);

/*======== void generate_sphere() ==========
Inputs:     struct matrix *points,
            double step,
            double x,
            double y,
            double z,
            double radius
Returns:
Adds the points of a sphere centered at (x, y, z) with a radius equal to the one
given to the matrix 'points' using the given step size (0 < step < 2).
====================================*/
void generate_sphere(struct matrix *points,
                     double step,
                     double x,
                     double y,
                     double z,
                     double radius);

/*======== void add_torus() ==========
Inputs:     struct matrix *points,
            double step,
            double x,
            double y,
            double z,
            double circle_radius,
            double torus_radius,
            drawing_mode draw_mode
Returns:
Adds the edges/polygons of a torus centered at (x, y, z) with a circle radius of
circle_radius and a torus radius of torus_radius using the given step size 
(0 < step < 2) to 'points'.
The draw_mode parameter determines the method used for drawing the prism. See
drawing_mode for more information.
====================================*/
void add_torus(struct matrix *points,
               double step,
               double x,
               double y,
               double z,
               double circle_radius,
               double torus_radius,
               drawing_mode draw_mode);

/*======== void generate_torus() ==========
Inputs:     struct matrix *points,
            double step,
            double x,
            double y,
            double z,
            double circle_radius,
            double torus_radius
Returns:
Adds the points of a torus centered at (x, y, z) with a circle radius of
circle_radius and a torus radius of torus_radius using the given step size 
(0 < step < 2) to 'points'.
====================================*/
void generate_torus(struct matrix *points,
                    double step,
                    double x,
                    double y,
                    double z,
                    double circle_radius,
                    double torus_radius);

/*======== void draw_line() ==========
Inputs:     screen s
            color c
            double x0
            double y0
            double z0
            double x1
            double y1
            double z1
            plotting_mode plot_mode
Returns:
Plots all the points necessary to draw line (x0, y0) - (x1, y1) onto
screen c using color c.
The plotting mode determines the coordinate system to be used when plotting points.
====================================*/
void draw_line(screen s, color c,
               double x0, double y0, double z0,
               double x1, double y1, double z1,
               plotting_mode plot_mode);

/*======== void draw_lines() ==========
Inputs:     screen s
            color c
            struct matrix * points
            plotting_mode plot_mode
Returns:
Iterates through 'points' 2 at a time and calls draw_line() to add that line
to the screen.
The plotting mode determines the coordinate system to be used when plotting points.
See plotting_mode for more information.
=====================================*/
void draw_lines(screen s, color c, struct matrix * points, plotting_mode plot_mode);

/*======== void draw_axes() ==========
Inputs:     screen s
            color c
Returns:
Plots all the points necessary to draw the x- and y-axes in the Cartesian
coordinate plane.
====================================*/
void draw_axes(screen s, color c);

/*======== void add_polygon() ==========
Inputs:     struct matrix *polygons
            double x0
            double y0
            double z0
            double x1
            double y1
            double z1
            double x2
            double y2
            double z2
Returns:
Adds the vertices (x0, y0, z0), (x1, y1, z1) and (x2, y2, z2) to the polygon
matrix. They define a single triangle surface.
======================================*/
void add_polygon(struct matrix *polygons,
                 double x0, double y0, double z0,
                 double x1, double y1, double z1,
                 double x2, double y2, double z2);


/*======== char is_visible() ==========
Inputs:     double **polygons
            int index
Returns:
Determines whether a polygon at the given index in the polygon matrix is
visible.  Uses the view vector and the normal vector to the polygon to calculate
the viewer's angle. Returns 1 if a polygon is visible and 0 if not.
======================================*/
char is_visible(double **polygons, int index);

/*======== void draw_polygons() ==========
Inputs:     screen s
            color c
            struct matrix *polygons
            plotting_mode plot_mode
Returns:
Iterates through the points in 'polygons' 3 points at a time, connecting them
to create triangles.
The plotting mode determines the coordinate system to be used when plotting points.
See plotting_mode for more information.
=========================================*/
void draw_polygons(screen s, color c, struct matrix *polygons, plotting_mode plot_mode);

/*======== void draw() ==========
Inputs:     screen s
            struct matrix *pts
            color c
Returns:
Draws the points in `pts` to the screen `s` using the color `c`.
This function automatically draws the axes if `global_plot_mode` is
PLOT_CARTESIAN, and it obeys `global_draw_mode` to either draw the points as
lines or polygons.
===============================*/
void draw(screen s, struct matrix *pts, color c);

/*======== void scanline_convert() ==========
Inputs:     screen s
            color c
            plotting_mode plot_mode
            double x0
            double y0
            double z0
            double x1
            double y1
            double z1
            double x2
            double y2
            double z2
Returns:
Performs scanline conversion on the triangle bounded by the three vertices 
(x0, y0, z0), (x1, y1, z1), and (x2, y2, z2).
============================================*/
void scanline_convert(screen s,
                      color c,
                      plotting_mode plot_mode,
                      double x0, double y0, double z0,
                      double x1, double y1, double z1,
                      double x2, double y2, double z2);
#endif
// vim: ts=4:et:sts:sw=4:sr
