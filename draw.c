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
    for (t = 0; t < 2 + step; t += step) {
        double rad = M_PI * (t + step);
        double new_x = cx + r * cos(rad);
        double new_y = cy + r * sin(rad);
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
    for (t = 0; t < 1 + step; t += step) {
        double t_squared = t * t;
        double t_cubed = t_squared * t;
        // Calculate new x and y coordinates using at^3 + bt^2 + ct + d
        double new_x = x_a * t_cubed + x_b * t_squared + x_c * t + x_d;
        double new_y = y_a * t_cubed + y_b * t_squared + y_c * t + y_d;
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
    for (t = 0; t < 1 + step; t += step) {
        double t_squared = t * t;
        double t_cubed = t_squared * t;
        // Calculate new x and y coordinates using at^3 + bt^2 + ct + d
        double new_x = x_a * t_cubed + x_b * t_squared + x_c * t + x_d;
        double new_y = y_a * t_cubed + y_b * t_squared + y_c * t + y_d;
        add_edge(points, old_x, old_y, 0, new_x, new_y, 0);
        old_x = new_x;
        old_y = new_y;
    }
    free_matrix(x_coeff);
    free_matrix(y_coeff);
}

void add_box(struct matrix *points,
             double x,
             double y,
             double z,
             double width,
             double height,
             double depth,
             drawing_mode draw_mode) {
    /* x, y, and z are the coordinates of the upper-left corner of the
     * front face of the rectangular prism, with width, height, and depth
     * corresponding to the x, y, and z coordinates respectively.
     */
    double x1 = x + width;
    double y1 = y - height;
    double z1 = z - depth;
    switch (draw_mode) {
        case DRAW_POLYGON:
            // Front side
            add_polygon(points, x, y, z, x, y1, z, x1, y1, z);
            add_polygon(points, x1, y1, z, x1, y, z, x, y, z);

            // Back side
            add_polygon(points, x, y, z1, x1, y, z1, x1, y1, z1);
            add_polygon(points, x1, y1, z1, x, y1, z1, x, y, z1);

            // Top side
            add_polygon(points, x, y, z1, x, y, z, x1, y, z);
            add_polygon(points, x1, y, z, x1, y, z1, x, y, z1);

            // Bottom side
            add_polygon(points, x, y1, z, x, y1, z1, x1, y1, z1);
            add_polygon(points, x1, y1, z1, x1, y1, z, x, y1, z);

            // Left side
            add_polygon(points, x, y, z, x, y, z1, x, y1, z1);
            add_polygon(points, x, y1, z1, x, y1, z, x, y, z);

            // Right side
            add_polygon(points, x1, y, z, x1, y1, z, x1, y1, z1);
            add_polygon(points, x1, y1, z1, x1, y, z1, x1, y, z);
            break;
        case DRAW_LINE:
            // Front side
            add_edge(points, x, y, z, x, y1, z);
            add_edge(points, x, y1, z, x1, y1, z);
            add_edge(points, x1, y1, z, x1, y, z);
            add_edge(points, x1, y, z, x, y, z);

            // Back side
            add_edge(points, x, y, z1, x, y1, z1);
            add_edge(points, x, y1, z1, x1, y1, z1);
            add_edge(points, x1, y1, z1, x1, y, z1);
            add_edge(points, x1, y, z1, x, y, z1);

            // Top side
            add_edge(points, x, y, z, x1, y, z);
            add_edge(points, x1, y, z, x1, y, z1);
            add_edge(points, x1, y, z1, x, y, z1);
            add_edge(points, x, y, z1, x, y, z);

            // Bottom side
            add_edge(points, x, y1, z, x1, y1, z);
            add_edge(points, x1, y1, z, x1, y1, z1);
            add_edge(points, x1, y1, z1, x, y1, z1);
            add_edge(points, x, y1, z1, x, y1, z);

            // Left side
            add_edge(points, x, y, z, x, y, z1);
            add_edge(points, x, y, z1, x, y1, z1);
            add_edge(points, x, y1, z1, x, y1, z);
            add_edge(points, x, y1, z, x, y, z);

            // Right side
            add_edge(points, x1, y, z, x1, y, z1);
            add_edge(points, x1, y, z1, x1, y1, z1);
            add_edge(points, x1, y1, z1, x1, y1, z);
            add_edge(points, x1, y1, z, x1, y, z);
            break;
    }
}

void add_sphere(struct matrix *points,
                double step,
                double x,
                double y,
                double z,
                double radius,
                drawing_mode draw_mode) {
    // Ensure that step size is greater than the minimum step size
    if (step < MIN_STEP_SIZE) {
        step = MIN_STEP_SIZE;
    }
    struct matrix *tmp = new_matrix(4, 1);
    generate_sphere(tmp, step, x, y, z, radius);
    // The points in the sphere are ordered as follows (assuming num_steps = 11):
    // P0  P1  P2  P3  P4  ... P10 P11
    // |A\B|                  
    // P12 P13 P14 P14 P15 ... P21 P22
    // P23 P24 P25 P26 P27 ... P32 P33
    // ...
    int latitude, longitude;
    double **m = tmp->m;
    int num_steps = round(1.0 / step) + 1;
    switch (draw_mode) {
        case DRAW_POLYGON: ; // Obligatory empty statement
            int num_pts = tmp->lastcol;
            for (latitude = 0; latitude < num_steps - 1; ++latitude) {
                int lat_start = num_steps * latitude;
                int next_lat_start = (lat_start + num_steps) % num_pts;
                for (longitude = 0; longitude < num_steps - 2; ++longitude) {
                    int index = lat_start + longitude;
                    int index_plus_one = index + 1;
                    int index_next_lat = next_lat_start + longitude;
                    int index_next_lat_plus_one = index_next_lat + 1;

                    /* DEBUG
                    print_debug("lat: %d, long: %d => %d", latitude, longitude, index);
                    print_debug("lastcol: %d", num_pts);
                    print_debug("index:%d", index);
                    print_debug("index_plus_one:%d", index_plus_one);
                    print_debug("index_next_lat:%d", index_next_lat);
                    print_debug("index_next_lat_plus_one:%d", index_next_lat_plus_one);
                    print_debug("index:\t\t\t\t(%lf,%lf,%lf)", m[0][index], m[1][index], m[2][index]);
                    print_debug("index_plus_one:\t\t\t(%lf,%lf,%lf)", m[0][index_plus_one], m[1][index_plus_one], m[2][index_plus_one]);
                    print_debug("index_next_lat:\t\t\t(%lf,%lf,%lf)", m[0][index_next_lat], m[1][index_next_lat], m[2][index_next_lat]);
                    print_debug("index_next_lat_plus_one:\t(%lf,%lf,%lf)", m[0][index_next_lat_plus_one], m[1][index_next_lat_plus_one], m[2][index_next_lat_plus_one]);
                    */

                    // Polygon-wise version
                    add_polygon(points,
                                m[0][index], m[1][index], m[2][index],
                                m[0][index_next_lat],
                                m[1][index_next_lat],
                                m[2][index_next_lat],
                                m[0][index_next_lat_plus_one],
                                m[1][index_next_lat_plus_one],
                                m[2][index_next_lat_plus_one]
                    );
                    add_polygon(points,
                                m[0][index_next_lat_plus_one],
                                m[1][index_next_lat_plus_one],
                                m[2][index_next_lat_plus_one],
                                m[0][index_plus_one],
                                m[1][index_plus_one],
                                m[2][index_plus_one],
                                m[0][index], m[1][index], m[2][index]
                    );
                }
            }
            break;
        case DRAW_LINE:
            // Line-wise version
            for (latitude = 0; latitude < num_steps; ++latitude) {
                for (longitude = 0; longitude < num_steps; ++longitude) {
                    int index = latitude * num_steps + longitude;
                    add_edge(points, m[0][index], m[1][index], m[2][index],
                                     m[0][index], m[1][index], m[2][index]);
                }
            }
            break;
    }
    free_matrix(tmp);
}

void generate_sphere(struct matrix *points,
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
    int t, s;
    int circle_steps, rotate_steps;
    circle_steps = rotate_steps = round(1.0 / step);
    for (s = 0; s < rotate_steps; ++s) {
        double phi_rad = M_PI * (2.0 * s / rotate_steps);
        double r_cos_phi = radius * cos(phi_rad);
        double r_sin_phi = radius * sin(phi_rad);
        for (t = 0; t <= circle_steps; ++t) {
            double theta_rad = M_PI * (1.0 * t / circle_steps);
            double sin_theta = sin(theta_rad);
            double cos_theta = cos(theta_rad);
            double curr_x = x + radius * cos_theta;
            double curr_y = y + sin_theta * r_cos_phi;
            double curr_z = z + sin_theta * r_sin_phi;
            add_point(points, curr_x, curr_y, curr_z);
        }
    }
}

void add_torus(struct matrix *points,
               double step,
               double x,
               double y,
               double z,
               double circle_radius,
               double torus_radius,
               drawing_mode draw_mode) {
    // Ensure that step size is greater than the minimum step size
    if (step < MIN_STEP_SIZE) {
        step = MIN_STEP_SIZE;
    }
    struct matrix *tmp = new_matrix(4, 1);
    generate_torus(tmp, step, x, y, z, circle_radius, torus_radius);
    // The points in the torus are ordered as follows (assuming num_steps = 10):
    // P0  P1  P2  P3  P4  ... P9
    // |A\B|                  
    // P10 P11 P12 P13 P14 ... P19
    // P20 P21 P22 P23 P24 ... P29
    // ...
    // P90 P91 P92 P93 P99 ... P99
    int latitude, longitude;
    double **m = tmp->m;
    int num_steps = round(1.0 / step);
    switch (draw_mode) {
        case DRAW_POLYGON: ; // Obligatory empty statement
            int num_pts = tmp->lastcol;
            for (latitude = 0; latitude < num_steps; ++latitude) {
                int lat_start = num_steps * latitude;
                int next_lat_start = (lat_start + num_steps) % num_pts;
                for (longitude = 0; longitude < num_steps; ++longitude) {
                    int index = lat_start + longitude;
                    // These checks ensure that the last point plotted is
                    // connected back to the correct starting point to close the
                    // loop/arc of the circle.
                    int index_plus_one = lat_start + ((longitude + 1) % num_steps);
                    int index_next_lat = (next_lat_start + longitude) % num_pts;
                    int index_next_lat_plus_one = (index_next_lat + 1) % num_pts;
                    add_polygon(points,
                                m[0][index], m[1][index], m[2][index],
                                m[0][index_next_lat],
                                m[1][index_next_lat],
                                m[2][index_next_lat],
                                m[0][index_next_lat_plus_one],
                                m[1][index_next_lat_plus_one],
                                m[2][index_next_lat_plus_one]
                    );
                    add_polygon(points,
                                m[0][index], m[1][index], m[2][index],
                                m[0][index_next_lat_plus_one],
                                m[1][index_next_lat_plus_one],
                                m[2][index_next_lat_plus_one],
                                m[0][index_plus_one],
                                m[1][index_plus_one],
                                m[2][index_plus_one]
                    );
                }
            }
            break;
        case DRAW_LINE:
            // Line-wise version
            for (latitude = 0; latitude < num_steps; ++latitude) {
                for (longitude = 0; longitude < num_steps; ++longitude) {
                    int index = latitude * num_steps + longitude;
                    add_edge(points, m[0][index], m[1][index], m[2][index],
                                     m[0][index], m[1][index], m[2][index]);
                }
            }
            break;
            break;
    }
    free_matrix(tmp);
}

void generate_torus(struct matrix *points,
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
          y-rotation       Circle points           Torus points
                           + translation
    where:
        Θ is the angle for generating the circle
        φ is the angle for generating the torus
        r is the radius of the circle
        R is the radius of the torus
    */
    int i, u;
    int torus_steps, circle_steps;
    torus_steps = circle_steps = round(1.0 / step);
    for (i = 0; i < circle_steps; ++i) {
        double theta_rad = M_PI * (2.0 * i / circle_steps);
        double circle_radius_cos_theta = circle_radius * cos(theta_rad);
        double circle_radius_sin_theta = circle_radius * sin(theta_rad);
        for (u = 0; u < torus_steps; ++u) {
            double phi_rad = M_PI * (2.0 * u / torus_steps);
            double curr_x = x + cos(phi_rad) * (circle_radius_cos_theta + torus_radius);
            double curr_y = y + circle_radius_sin_theta;
            double curr_z = z - sin(phi_rad) * (circle_radius_cos_theta + torus_radius);
            add_point(points, curr_x, curr_y, curr_z);
        }
    }
}

void draw_line(screen s, color c, double x0, double y0, double x1, double y1,
               plotting_mode plot_mode) {
    // Ensure that x values are increasing (or equal), for simplification
    if (x0 > x1) {
        // Swap points if necessary
        double tmp;
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

void draw_lines(screen s, color c, struct matrix *edges,
                plotting_mode plot_mode) {
    if (edges->lastcol < 2) {
        print_error("Edge matrix has a length of less than 2.");
        return;
    }
    int i;
    for (i = 0; i < edges->lastcol - 1; i+=2) {
        draw_line(s, c, edges->m[0][i], edges->m[1][i],
                  edges->m[0][i+1], edges->m[1][i+1], plot_mode);
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

void add_polygon(struct matrix *polygons,
                 double x0, double y0, double z0,
                 double x1, double y1, double z1,
                 double x2, double y2, double z2) {
    add_point(polygons, x0, y0, z0);
    add_point(polygons, x1, y1, z1);
    add_point(polygons, x2, y2, z2);
}

char is_visible(double **polygons, int index) {
    // cos(Θ) = (N • V) / (|N| * |V|)
    // If cos(Θ) < 0, then 90° < Θ < 270°
    // Supplementary angle: 0° < Φ < 90°, so the surface is visible
    double *normal = cross_prod(polygons[0][index+1] - polygons[0][index],
                                polygons[1][index+1] - polygons[1][index],
                                polygons[2][index+1] - polygons[2][index],
                                polygons[0][index+2] - polygons[0][index],
                                polygons[1][index+2] - polygons[1][index],
                                polygons[2][index+2] - polygons[2][index]);
    double view_vector[3] = {0.0, 0.0, -1.0};
    char visible = (dot_prod(normal[0], normal[1], normal[2],
                             view_vector[0], view_vector[1], view_vector[2]
                           ) < 0) ? 1 : 0;
    free(normal);
    return visible;
}

void draw_polygons(screen s, color c, struct matrix *polygons,
                   plotting_mode plot_mode) {
    if (polygons->lastcol < 3) {
        print_error("Polygon matrix has a length of less than 3.");
        return;
    }
    double **m = polygons->m;
    int i;
    for (i = 0; i < polygons->lastcol - 2; i+=3) {
        if (is_visible(m, i)) {
            draw_line(s, c, m[0][i], m[1][i],
                            m[0][i+1], m[1][i+1],
                            plot_mode);
            draw_line(s, c, m[0][i+1], m[1][i+1],
                            m[0][i+2], m[1][i+2],
                            plot_mode);
            draw_line(s, c, m[0][i+2], m[1][i+2],
                            m[0][i], m[1][i],
                            plot_mode);
            /* DEBUG
            print_debug("drawing (%lf,%lf,%lf)\n"
                "\tto      (%lf, %lf, %lf)\n"
                "\tand     (%lf, %lf, %lf)",
                m[0][i], m[1][i], m[2][i],
                m[0][i+1], m[1][i+1], m[2][i+1],
                m[0][i+2], m[1][i+2], m[2][i+2]);
            if (i % 20 == 0) {
                display(s);
            }
            */
        }
    }
}

// vim: ts=4:et:sts:sw=4:sr
