#include "draw.h"

void add_point(struct matrix * points, double x, double y, double z) {
    if (points->lastcol >= points->cols - 1) {
        grow_matrix(points, points->cols * 2);
    }
    double **m = points->m;
    int lastcol = points->lastcol;
    m[0][lastcol] = x;
    m[1][lastcol] = y;
    m[2][lastcol] = z;
    m[3][lastcol] = 1;
    ++points->lastcol;
}

void add_edge(struct matrix * points,
              double x0, double y0, double z0,
              double x1, double y1, double z1) {
    add_point(points, x0, y0, z0);
    add_point(points, x1, y1, z1);
}

void add_circle(struct matrix *points,
                double cx,
                double cy,
                double cz,
                double r,
                double step) {
    double t;
    double old_x = cx + r; // r*cos(0) is r
    double old_y = cy + 0; // r*sin(0) is 0
    double new_x;
    double new_y;
    for (t = 0; t < 2 + step; t += step) {
        double rad = M_PI * (t + step);
        new_x = cx + r * cos(rad);
        new_y = cy + r * sin(rad);
        add_edge(points, old_x, old_y, cz, new_x, new_y, cz);
        old_x = new_x;
        old_y = new_y;
    }
}

void add_curve(struct matrix *points,
               double step,
               curve_type type,
               double x0, double y0,
               double x1, double y1,
               double x2, double y2,
               double x3, double y3) {
    switch (type) {
        case HERMITE_CURVE:
            add_hermite_curve(points, step,
                            x0, y0,
                            x2, y2,
                            x1 - x0, y1 - y0,
                            x3 - x2, y3 - y2);
            break;
        case BEZIER_CURVE:
            add_bezier_curve(points, step,
                                x0, y0,
                                x1, y1,
                                x2, y2,
                                x3, y3);
            break;
    }
}

void add_hermite_curve(struct matrix *points,
                       double step,
                       double x0, double y0,
                       double x1, double y1,
                       double dx0, double dy0,
                       double dx1, double dy1) {
    /*
    |  2 -2  1  1 || P0 |   | a |
    | -3  3 -2 -1 || P1 | = | b |
    |  0  0  1  0 || R0 |   | c |
    |  1  0  0  0 || R1 |   | d |
    */
    struct matrix *coeff = make_hermite_coefficients();
    struct matrix *x = new_matrix(4, 1);
    x->m[0][0] = x0;
    x->m[1][0] = x1;
    x->m[2][0] = dx0;
    x->m[3][0] = dx1;
    x->lastcol = 1;
    // Calculate a, b, c, and d in at^3 + bt^2 + ct + d
    // used to determine the x coordinate of the next point
    struct matrix *x_coeff = matrix_mult(coeff, x);
    struct matrix *y = new_matrix(4, 1);
    y->m[0][0] = y0;
    y->m[1][0] = y1;
    y->m[2][0] = dy0;
    y->m[3][0] = dy1;
    y->lastcol = 1;
    // Calculate a, b, c, and d in at^3 + bt^2 + ct + d
    // used to determine the y coordinate of the next point
    struct matrix *y_coeff = matrix_mult(coeff, y);
    free_matrix(coeff);
    free_matrix(x);
    free_matrix(y);
    double x_a = x_coeff->m[0][0];
    double x_b = x_coeff->m[1][0];
    double x_c = x_coeff->m[2][0];
    double x_d = x_coeff->m[3][0];
    double y_a = y_coeff->m[0][0];
    double y_b = y_coeff->m[1][0];
    double y_c = y_coeff->m[2][0];
    double y_d = y_coeff->m[3][0];
    double t;
    double old_x = x0;
    double old_y = y0;
    double new_x;
    double new_y;
    for (t = 0; t < 1 + step; t += step) {
        double t_squared = t * t;
        double t_cubed = t_squared * t;
        // Calculate new x and y coordinates using at^3 + bt^2 + ct + d
        new_x = x_a * t_cubed + x_b * t_squared + x_c * t + x_d;
        new_y = y_a * t_cubed + y_b * t_squared + y_c * t + y_d;
        add_edge(points, old_x, old_y, 0, new_x, new_y, 0);
        old_x = new_x;
        old_y = new_y;
    }
    free_matrix(x_coeff);
    free_matrix(y_coeff);
}

void add_bezier_curve(struct matrix *points,
                      double step,
                      double x0, double y0,
                      double x1, double y1,
                      double x2, double y2,
                      double x3, double y3) {
    /*
    | -1  3 -3 1 |   | P0 |   | a |
    |  3 -6  3 0 | * | P1 | = | b |
    | -3  3  0 0 |   | P2 |   | c |
    |  1  0  0 0 |   | P3 |   | d |
    */
    struct matrix *coeff = make_bezier_coefficients();
    struct matrix *x = new_matrix(4, 1);
    x->m[0][0] = x0;
    x->m[1][0] = x1;
    x->m[2][0] = x2;
    x->m[3][0] = x3;
    x->lastcol = 1;
    // Calculate a, b, c, and d in at^3 + bt^2 + ct + d
    // used to determine the x coordinate of the next point
    struct matrix *x_coeff = matrix_mult(coeff, x);
    struct matrix *y = new_matrix(4, 1);
    y->m[0][0] = y0;
    y->m[1][0] = y1;
    y->m[2][0] = y2;
    y->m[3][0] = y3;
    y->lastcol = 1;
    // Calculate a, b, c, and d in at^3 + bt^2 + ct + d
    // used to determine the y coordinate of the next point
    struct matrix *y_coeff = matrix_mult(coeff, y);
    free_matrix(coeff);
    free_matrix(x);
    free_matrix(y);
    double x_a = x_coeff->m[0][0];
    double x_b = x_coeff->m[1][0];
    double x_c = x_coeff->m[2][0];
    double x_d = x_coeff->m[3][0];
    double y_a = y_coeff->m[0][0];
    double y_b = y_coeff->m[1][0];
    double y_c = y_coeff->m[2][0];
    double y_d = y_coeff->m[3][0];
    double t;
    double old_x = x0;
    double old_y = y0;
    double new_x;
    double new_y;
    for (t = 0; t < 1 + step; t += step) {
        double t_squared = t * t;
        double t_cubed = t_squared * t;
        // Calculate new x and y coordinates using at^3 + bt^2 + ct + d
        new_x = x_a * t_cubed + x_b * t_squared + x_c * t + x_d;
        new_y = y_a * t_cubed + y_b * t_squared + y_c * t + y_d;
        add_edge(points, old_x, old_y, 0, new_x, new_y, 0);
        old_x = new_x;
        old_y = new_y;
    }
    free_matrix(x_coeff);
    free_matrix(y_coeff);
}

void add_prism(struct matrix *points,
               double x,
               double y,
               double z,
               double width,
               double height,
               double depth) {
    /* x, y, and z are the coordinates of the upper-left corner of the
     * front face of the rectangular prism, with width, height, and depth
     * corresponding to the x, y, and z coordinates respectively.
     */
    add_edge(points, x, y, z, x, y, z);
    add_edge(points, x + width, y, z, x + width, y, z);
    add_edge(points, x + width, y - height, z, x + width, y - height, z);
    add_edge(points, x, y - height, z, x, y - height, z);
    add_edge(points, x, y, z - depth, x, y, z - depth);
    add_edge(points, x + width, y, z - depth, x + width, y, z - depth);
    add_edge(points, x + width, y - height, z - depth, x + width, y - height, z - depth);
    add_edge(points, x, y - height, z - depth, x, y - height, z - depth);
}

void add_sphere(struct matrix *points,
                double step,
                double x,
                double y,
                double z,
                double radius) {
    /*
    | 1     0    0     1 || rcos(Θ) |   |    rcos(Θ)    |
    | 1 cos(φ) -sin(φ) 0 || rsin(Θ) |   | rsin(Θ)cos(φ) |
    | 0 sin(φ)  cos(φ) 0 ||    0    | = | rsin(Θ)sin(φ) |
    | 0    0     0     1 ||    1    |   |       1       |
          x-rotation     Circle points    Sphere points
    where:
        Θ is the angle for generating the circle
        φ is the angle for generating the sphere
        r is the radius of the sphere
    */
    // Ensure that step size is greater than the minimum step size
    if (step < MIN_STEP_SIZE) {
        step = MIN_STEP_SIZE;
    }
    double old_x = x + radius; // rcos(0) = r
    double old_y = y;          // rsin(0)cos(0) = 0
    double old_z = z;          // rsin(0)sin(0) = 0
    double t, s;
    for (t = 0; t < 2 + step; t += step) {
        double theta_rad = M_PI * (t + step);
        double r_cos_theta = radius * cos(theta_rad);
        double r_sin_theta = radius * sin(theta_rad);
        for (s = 0; s < 1 + step; s += step) {
            double phi_rad = M_PI * (s + step);
            double new_x = x + r_cos_theta;
            double new_y = y + r_sin_theta * cos(phi_rad);
            double new_z = z + r_sin_theta * sin(phi_rad);
            //add_edge(points, old_x, old_y, old_z, new_x, new_y, new_z);
            add_edge(points, old_x, old_y, old_z, old_x, old_y, old_z);
            old_x = new_x;
            old_y = new_y;
            old_z = new_z;
        }
    }
}

void add_torus(struct matrix *points,
               double step,
               double x,
               double y,
               double z,
               double circle_radius,
               double torus_radius) {
    /*
    | cos(φ) 0 -sin(φ) 0 || rcos(Θ) + R |   | cos(φ) * (rcos(Θ) + R) |
    | 0    1     0     0 ||   rsin(Θ)   |   |         rsin(Θ)        |
    | sin(φ) 0 cos(φ)  0 ||      0      | = |  -sin(φ)(rcos(Θ) + R)  |
    | 0    0     0     1 ||      1      |   |            1           |
          y-rotation       Circle points           Sphere points
                           + translation
    where:
        Θ is the angle for generating the circle
        φ is the angle for generating the torus
        r is the radius of the circle
        R is the radius of the torus
    */
    // Ensure that step size is greater than the minimum step size
    if (step < MIN_STEP_SIZE) {
        step = MIN_STEP_SIZE;
    }
    double old_x = x + circle_radius + torus_radius; // cos(0) * (rcos(0) + R) = r + R
    double old_y = y;                                // rsin(0) = 0
    double old_z = z;                                // -sin(0)(rcos(0) + R) = 0
    double t, s;
    for (t = 0; t < 2 + step; t += step) {
        double theta_rad = M_PI * (t + step);
        double circle_radius_cos_theta = circle_radius * cos(theta_rad);
        double circle_radius_sin_theta = circle_radius * sin(theta_rad);
        for (s = 0; s < 2 + step; s += step) {
            double phi_rad = M_PI * (s + step);
            double new_x = x + cos(phi_rad) * (circle_radius_cos_theta + torus_radius);
            double new_y = y + circle_radius_sin_theta;
            double new_z = z - sin(phi_rad) * (circle_radius_cos_theta + torus_radius);
            //add_edge(points, old_x, old_y, old_z, new_x, new_y, new_z);
            add_edge(points, old_x, old_y, old_z, old_x, old_y, old_z);
            old_x = new_x;
            old_y = new_y;
            old_z = new_z;
        }
    }
}

void draw_line(screen s, color c, double x0, double y0, double x1, double y1,
               plotting_mode plot_mode) {
    // Ensure that x values are increasing (or equal), for simplification
    double tmp;
    if (x0 > x1) {
        // Swap points if necessary
        tmp = x0;
        x0 = x1;
        x1 = tmp;
        tmp = y0;
        y0 = y1;
        y1 = tmp;
    }
    double a = 2 * (y1 - y0);
    double b = -2 * (x1 - x0);
    // 1st and 5th octants of the 2D plane
    if (a >= 0 && a <= (-1*b)) {
        // Shortcut for the first round of calculations of Ax + By + C
        // using the midpoint of the next two possible coordinates:
        // d = f(x0+1, y0+1/2)
        //   = A(x0+1) + B(y+1/2) + C
        //   = Ax0 + By0 + C + A + B/2
        //   = A + B/2
        double d = a + b / 2;
        double x = x0;
        double y = y0;
        while (x <= x1) {
            switch (plot_mode) {
                case PLOT_CARTESIAN:
                    plot_cartesian(s, c, x, y);
                    break;
                case PLOT_ABSOLUTE:
                    plot_absolute(s, c, x, y);
                    break;
            }
            // If Ax + By + C > 0, then the midpoint of the next two
            // possible coordinates is below the line,
            // so we have to draw the pixel in the upper coordinate
            // We increase d according to the rule:
            // if   x → x + 1
            //      y → y
            // then d = d + A
            // if   x → x + 1
            //      y → y + 1
            // then d = d + A + B
            if (d > 0) {
                ++y;
                d += b;
            }
            ++x;
            d += a;
        }
    }
    // 2nd and 6th octants of the 2D plane
    else if (a >= 0 && a > (-1*b)) {
        // Shortcut for the first round of calculations of Ax + By + C
        // using the midpoint of the next two possible coordinates:
        // d = f(x0+1/2, y0+1)
        //   = A(x0+1/2) + B(y+1) + C
        //   = Ax0 + By0 + C + A/2 + B
        //   = A/2 + B
        double d = a / 2 + b;
        double x = x0;
        double y = y0;
        while (y <= y1) {
            switch (plot_mode) {
                case PLOT_CARTESIAN:
                    plot_cartesian(s, c, x, y);
                    break;
                case PLOT_ABSOLUTE:
                    plot_absolute(s, c, x, y);
                    break;
            }
            // If Ax + By + C < 0, then the midpoint of the next two
            // possible coordinates is above the line,
            // so we have to draw the pixel in the righter coordinate
            // We increase d according to the rule:
            // if   y → y + 1
            //      x → x
            // then d = d + B
            // if   y → y + 1
            //      x → x + 1
            // then d = d + A + B
            if (d < 0) {
                ++x;
                d += a;
            }
            ++y;
            d += b;
        }
    }
    // 3rd and 7th octants of the 2D plane
    else if (a < 0 && a <= b) {
        // Shortcut for the first round of calculations of Ax + By + C
        // using the midpoint of the next two possible coordinates:
        // d = f(x0+1/2, y0-1)
        //   = A(x0+1/2) + B(y0-1) + C
        //   = Ax0 + By0 + C + A/2 - B
        //   = A/2 - B
        double d = a / 2 - b;
        double x = x0;
        double y = y0;
        while (y >= y1) {
            switch (plot_mode) {
                case PLOT_CARTESIAN:
                    plot_cartesian(s, c, x, y);
                    break;
                case PLOT_ABSOLUTE:
                    plot_absolute(s, c, x, y);
                    break;
            }
            // If Ax + By + C > 0, then the midpoint of the next two
            // possible coordinates is below the line,
            // so we have to draw the pixel in the righter coordinate
            // We increase d according to the rule:
            // if   y → y - 1
            //      x → x
            // then d = d - B
            // if   y → y - 1
            //      x → x + 1
            // then d = d + A - B
            if (d > 0) {
                ++x;
                d += a;
            }
            --y;
            d -= b;
        }
    }
    // 4th and 8th octants of the 2D plane
    else if (a < 0 && a > b) {
        // Shortcut for the first round of calculations of Ax + By + C
        // using the midpoint of the next two possible coordinates:
        // d = f(x0+1, y0-1/2)
        //   = A(x0+1) + B(y0-1/2) + C
        //   = Ax0 + By0 + C + A - B/2
        //   = A - B/2
        double d = a - b / 2;
        double x = x0;
        double y = y0;
        while (x <= x1) {
            switch (plot_mode) {
                case PLOT_CARTESIAN:
                    plot_cartesian(s, c, x, y);
                    break;
                case PLOT_ABSOLUTE:
                    plot_absolute(s, c, x, y);
                    break;
            }
            // If Ax + By + C < 0, then the midpoint of the next two
            // possible coordinates is above the line,
            // so we have to draw the pixel in the lower coordinate
            // We increase d according to the rule:
            // if   x → x + 1
            //      y → y
            // then d = d + A
            // if   x → x + 1
            //      y → y - 1
            // then d = d + A - B
            if (d < 0) {
                --y;
                d -= b;
            }
            ++x;
            d += a;
        }
    }
}

void draw_lines(screen s, color c, struct matrix * points, plotting_mode plot_mode) {
    int i;
    for (i = 0; i < points->lastcol - 1; i+=2) {
        draw_line(s, c, points->m[0][i], points->m[1][i],
                  points->m[0][i+1], points->m[1][i+1], plot_mode);
    }
}

void draw_axes(screen s, color c) {
    int i;
    for (i = -1 * XRES_CARTESIAN; i < XRES_CARTESIAN; ++i) {
        plot_cartesian(s, c, i, 0);
    }
    for (i = -1 * YRES_CARTESIAN; i < YRES_CARTESIAN; ++i) {
        plot_cartesian(s, c, 0, i);
    }
}

// vim: ts=4:et:sts:sw=4:sr
