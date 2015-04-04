#include "parser.h"

FILE *global_file = NULL;
struct matrix *global_trans_mat = NULL;
struct matrix *global_pts = NULL;
screen global_s = NULL;
int line_no = 0;

void parse_input(char *cmd, char *line_buf, char error_is_fatal) {
    if (strcmp(cmd, LINE_CMD) == 0) {
        #ifdef DEBUG
        print_debug("Got line command");
        #endif
        double x1, x2, y1, y2, z1, z2;
        int retVal = sscanf(line_buf, "%*s %lf %lf %lf %lf %lf %lf",
                            &x1, &y1, &z1, &x2, &y2, &z2);
        if (retVal != 6) {
            print_error("Invalid arguments for line command: \"%s\"", line_buf);
            if (error_is_fatal) {
                exit(EXIT_FAILURE);
            }
            return;
        }
        add_edge(global_pts, x1, y1, z1, x2, y2, z2);
    }
    else if (strcmp(cmd, IDENTITY_CMD) == 0) {
        #ifdef DEBUG
        print_debug("Got identity command");
        #endif
        ident(global_trans_mat);
    }
    else if (strcmp(cmd, SCALE_CMD) == 0) {
        #ifdef DEBUG
        print_debug("Got scale command");
        #endif
        double sx, sy, sz;
        int retVal = sscanf(line_buf, "%*s %lf %lf %lf", &sx, &sy, &sz);
        if (retVal != 3) {
            print_error("Invalid arguments for scale command: %s", line_buf);
            if (error_is_fatal) {
                exit(EXIT_FAILURE);
            }
            return;
        }
        struct matrix *scale = make_scale(sx, sy, sz);
        struct matrix *tmp = apply_transform(scale, global_trans_mat);
        free_matrix(global_trans_mat);
        free_matrix(scale);
        global_trans_mat = tmp;
    }
    else if (strcmp(cmd, TRANSLATE_CMD) == 0) {
        #ifdef DEBUG
        print_debug("Got translate command");
        #endif
        double tx, ty, tz;
        int retVal = sscanf(line_buf, "%*s %lf %lf %lf", &tx, &ty, &tz);
        if (retVal != 3) {
            print_error("Invalid arguments for translate command: \"%s\"", line_buf);
            if (error_is_fatal) {
                exit(EXIT_FAILURE);
            }
            return;
        }
        struct matrix *translate = make_translate(tx, ty, tz);
        struct matrix *tmp = apply_transform(translate, global_trans_mat);
        free_matrix(global_trans_mat);
        free_matrix(translate);
        global_trans_mat = tmp;
    }
    else if (strcmp(cmd, XROT_CMD) == 0) {
        #ifdef DEBUG
        print_debug("Got x rotate command");
        #endif
        double theta;
        int retVal = sscanf(line_buf, "%*s %lf", &theta);
        if (retVal != 1) {
            print_error("Invalid arguments for x rotate command: %s", line_buf);
            if (error_is_fatal) {
                exit(EXIT_FAILURE);
            }
            return;
        }
        struct matrix *rotX = make_rotX(theta);
        struct matrix *tmp = apply_transform(rotX, global_trans_mat);
        free_matrix(global_trans_mat);
        free_matrix(rotX);
        global_trans_mat = tmp;
    }
    else if (strcmp(cmd, YROT_CMD) == 0) {
        #ifdef DEBUG
        print_debug("Got y rotate command");
        #endif
        double theta;
        int retVal = sscanf(line_buf, "%*s %lf", &theta);
        if (retVal != 1) {
            print_error("Invalid arguments for y rotate command: \"%s\"", line_buf);
            if (error_is_fatal) {
                exit(EXIT_FAILURE);
            }
            return;
        }
        struct matrix *rotY = make_rotY(theta);
        struct matrix *tmp = apply_transform(rotY, global_trans_mat);
        free_matrix(global_trans_mat);
        free_matrix(rotY);
        global_trans_mat = tmp;
    }
    else if (strcmp(cmd, ZROT_CMD) == 0) {
        #ifdef DEBUG
        print_debug("Got z rotate command");
        #endif
        double theta;
        int retVal = sscanf(line_buf, "%*s %lf", &theta);
        if (retVal != 1) {
            print_error("Invalid arguments for z rotate command: \"%s\"", line_buf);
            if (error_is_fatal) {
                exit(EXIT_FAILURE);
            }
            return;
        }
        struct matrix *rotZ = make_rotZ(theta);
        struct matrix *tmp = apply_transform(rotZ, global_trans_mat);
        free_matrix(global_trans_mat);
        free_matrix(rotZ);
        global_trans_mat = tmp;
    }
    else if (strcmp(cmd, APPLY_CMD) == 0) {
        #ifdef DEBUG
        print_debug("Got apply command");
        #endif
        struct matrix *tmp = apply_transform(global_trans_mat, global_pts);
        free_matrix(global_pts);
        global_pts = tmp;
    }
    else if (strcmp(cmd, CIRCLE_CMD) == 0) {
        #ifdef DEBUG
        print_debug("Got circle command");
        #endif
        double cx, cy, cz, r;
        int retVal = sscanf(line_buf, "%*s %lf %lf %lf %lf", &cx, &cy, &cz, &r);
        if (retVal != 4) {
            print_error("Invalid arguments for circle command: \"%s\"", line_buf);
            if (error_is_fatal) {
                exit(EXIT_FAILURE);
            }
            return;
        }
        add_circle(global_pts, cx, cy, cz, r, STEP_SIZE);
    }
    else if (strcmp(cmd, HERMITE_CURVE_CMD) == 0) {
        #ifdef DEBUG
        print_debug("Got Hermite curve command");
        #endif
        double x0, y0, x1, y1, x2, y2, x3, y3;
        int retVal = sscanf(line_buf, "%*s %lf %lf %lf %lf %lf %lf %lf %lf",
                            &x0, &y0, &x1, &y1, &x2, &y2, &x3, &y3);
        if (retVal != 8) {
            print_error("Invalid arguments for Hermite curve command: \"%s\"", line_buf);
            if (error_is_fatal) {
                exit(EXIT_FAILURE);
            }
            return;
        }
        add_curve(global_pts, STEP_SIZE, HERMITE_CURVE,
                  x0, y0, x1, y1, x2, y2, x3, y3);
    }
    else if (strcmp(cmd, BEZIER_CURVE_CMD) == 0) {
        #ifdef DEBUG
        print_debug("Got Bezier curve command");
        #endif
        double x0, y0, x1, y1, x2, y2, x3, y3;
        int retVal = sscanf(line_buf, "%*s %lf %lf %lf %lf %lf %lf %lf %lf",
                            &x0, &y0, &x1, &y1, &x2, &y2, &x3, &y3);
        if (retVal != 8) {
            print_error("Invalid arguments for Bezier curve command: \"%s\"", line_buf);
            if (error_is_fatal) {
                exit(EXIT_FAILURE);
            }
            return;
        }
        add_curve(global_pts, STEP_SIZE, BEZIER_CURVE,
                  x0, y0, x1, y1, x2, y2, x3, y3);
    }
    else if (strcmp(cmd, VIEW_CMD) == 0) {
        #ifdef DEBUG
        print_debug("Got view command");
        #endif
        // Draw image
        // First, draw the axes if plotting in Cartesian plane
        draw_axes_if_cartesian_mode();
        draw();
        // Display image
        display(global_s);
    }
    else if (strcmp(cmd, SAVE_CMD) == 0) {
        #ifdef DEBUG
        print_debug("Got save command");
        #endif
        // Draw image
        // First, draw the axes if plotting in Cartesian plane
        draw_axes_if_cartesian_mode();
        draw();
        // Save image to file with given filename
        save(line_buf, error_is_fatal);
    }
    else if (strcmp(cmd, VIEW_AND_SAVE_CMD) == 0) {
        #ifdef DEBUG
        print_debug("Got view and save command");
        #endif
        // Draw image
        // First, draw the axes if plotting in Cartesian plane
        draw_axes_if_cartesian_mode();
        draw();
        // Display image
        display(global_s);
        // Save image to file with given filename
        save(line_buf, error_is_fatal);
    }
    else if (strcmp(cmd, CLEAR_EDGE_MATRIX_CMD) == 0) {
        // Clear the point matrix
        #ifdef DEBUG
        print_debug("Got clear point matrix command");
        #endif
        free_matrix(global_pts);
        global_pts = new_matrix(4, 1);
    }
    else if (strcmp(cmd, PRISM_CMD) == 0) {
        // Add rectangular prism to edge matrix
        #ifdef DEBUG
        print_debug("Got draw prism command");
        #endif
        double x, y, z, width, height, depth;
        int retVal = sscanf(line_buf, "%*s %lf %lf %lf %lf %lf %lf",
                            &x, &y, &z, &width, &height, &depth);
        if (retVal != 6) {
            print_error("Invalid arguments for draw prism command: \"%s\"", line_buf);
            if (error_is_fatal) {
                exit(EXIT_FAILURE);
            }
            return;
        }
        add_prism(global_pts, x, y, z, width, height, depth, global_draw_mode);
    }
    else if (strcmp(cmd, SPHERE_CMD) == 0) {
        // Add sphere to edge matrix
        #ifdef DEBUG
        print_debug("Got draw sphere command");
        #endif
        double x, y, z, radius, step;
        int retVal = sscanf(line_buf, "%*s %lf %lf %lf %lf %lf",
                            &x, &y, &z, &radius, &step);
        if (retVal != 5) {
            print_error("Invalid arguments for draw sphere command: \"%s\"", line_buf);
            if (error_is_fatal) {
                exit(EXIT_FAILURE);
            }
            return;
        }
        add_sphere(global_pts, step, x, y, z, radius, global_draw_mode);
    }
    else if (strcmp(cmd, TORUS_CMD) == 0) {
        // Add torus to edge matrix
        #ifdef DEBUG
        print_debug("Got draw torus command");
        #endif
        double x, y, z, circleRadius, torusRadius, step;
        int retVal = sscanf(line_buf, "%*s %lf %lf %lf %lf %lf %lf",
                            &x, &y, &z, &circleRadius, &torusRadius, &step);
        if (retVal != 6) {
            print_error("Invalid arguments for draw torus command: \"%s\"", line_buf);
            if (error_is_fatal) {
                exit(EXIT_FAILURE);
            }
            return;
        }
        add_torus(global_pts, step, x, y, z,
                  circleRadius, torusRadius, global_draw_mode);
    }
    else if (strcmp(cmd, PRINT_TRANSFORMATION_MATRIX_CMD) == 0) {
        #ifdef DEBUG
        print_debug("Got print transformation matrix command");
        #endif
        print_matrix(global_trans_mat);
    }
    else if (strcmp(cmd, PRINT_POINT_MATRIX_CMD) == 0) {
        #ifdef DEBUG
        print_debug("Got print point matrix command");
        #endif
        print_matrix(global_pts);
    }
    else if (strcmp(cmd, PLOT_MODE_CMD) == 0) {
        #ifdef DEBUG
        print_debug("Got plot mode command");
        #endif
        char mode[10];
        int retVal = sscanf(line_buf, "%*s %9s", mode);
        if (strcmp(mode, CARTESIAN) == 0) {
            global_plot_mode = PLOT_CARTESIAN;
        }
        else if (strcmp(mode, ABSOLUTE) == 0) {
            global_plot_mode = PLOT_ABSOLUTE;
        }
        else {
            print_error("Invalid arguments for plot mode command: \"%s\"", line_buf);
            if (error_is_fatal) {
                exit(EXIT_FAILURE);
            }
            return;
        }
    }
    else if (strcmp(cmd, DRAW_MODE_CMD) == 0) {
        #ifdef DEBUG
        print_debug("Got draw mode command");
        #endif
        char mode[10];
        int retVal = sscanf(line_buf, "%*s %9s", mode);
        if (strcmp(mode, LINE) == 0) {
            global_draw_mode = DRAW_LINE;
        }
        else if (strcmp(mode, POLYGON) == 0) {
            global_draw_mode = DRAW_POLYGON;
        }
        else {
            print_error("Invalid arguments for plot mode command: \"%s\"", line_buf);
            if (error_is_fatal) {
                exit(EXIT_FAILURE);
            }
            return;
        }
    }
    else if (strcmp(cmd, QUIT_CMD) == 0) {
        #ifdef DEBUG
        print_debug("Got quit command");
        #endif
        free_matrix(global_pts);
        free_matrix(global_trans_mat);
        free_screen(global_s);
        fclose(global_file);
        exit(EXIT_SUCCESS);
    }
    else {
        print_error("Invalid command on line %d: \"%s\"", line_no, cmd);
    }
}

void parse_file(char error_is_fatal) {
    if (global_file == NULL
        || global_trans_mat == NULL
        || global_pts == NULL
        || global_s == NULL) {
        print_error("parse_file(): Global variables are not synchronized."
                    " Have you called synchronize_variables()?");
        return;
    }
    char line_buf[LINE_BUF_SIZE];
    line_no = 0;
    while (fgets(line_buf, sizeof(line_buf), global_file) != NULL) {
        char cmd[CMD_BUF_SIZE];
        int retVal = sscanf(line_buf, "%10s", cmd);
        // Skip empty input
        if (retVal < 0) {
            continue;
        }
        if (remove_trailing_newline(line_buf) == 0) {
            ++line_no;
        }
        #ifdef DEBUG
        print_debug("Line %d: \"%s\"", line_no, line_buf);
        #endif
        parse_input(cmd, line_buf, error_is_fatal);
    }
    fclose(global_file);
}

void synchronize_variables(FILE *file,
                           struct matrix *trans_mat,
                           struct matrix *pts,
                           screen s) {
    if (file != NULL) {
        global_file = file;
    }
    if (trans_mat != NULL) {
        if (trans_mat->rows != 4 && trans_mat->cols != 4) {
            print_error("Invalid dimensions for transformation matrix. "
                        "Given: (%d x %d); Required: (4 x 4)",
                        trans_mat->rows,
                        trans_mat->cols);
            exit(EXIT_FAILURE);
        }
        global_trans_mat = trans_mat;
        ident(global_trans_mat);
    }
    if (pts != NULL) {
        global_pts = pts;
    }
    if (s != NULL) {
        global_s = s;
    }
}

void free_variables() {
    #ifdef DEBUG
    print_debug("Freeing all alloc'd memory");
    #endif
    free_matrix(global_trans_mat);
    free_matrix(global_pts);
    free_screen(global_s);
    global_file = NULL;
    global_trans_mat = NULL;
    global_pts = NULL;
    global_s = NULL;
}

void draw_axes_if_cartesian_mode() {
    // Draw the axes if plotting in Cartesian plane
    color c;
    if (global_plot_mode == PLOT_CARTESIAN) {
        clear_screen(global_s);
        c.red = 255;
        c.blue = 0;
        c.green = 0;
        draw_axes(global_s, c);
    }
}

void draw() {
    // Draw the points matrix to the screen using the current drawing mode
    color c;
    c.red = 30;
    c.blue = 100;
    c.green = 155;
    switch (global_draw_mode) {
        case DRAW_LINE:
            draw_lines(global_s, c, global_pts, global_plot_mode);
            break;
        case DRAW_POLYGON:
            draw_polygons(global_s, c, global_pts, global_plot_mode);
            break;
    }
}

int save(char *line_buf, int error_is_fatal) {
    int retVal;
    char filename[101];
    retVal = sscanf(line_buf, "%*s %100s", filename);
    if (retVal == 1) {
        save_ppm(global_s, filename);
        save_extension(global_s, filename);
    }
    else {
        print_error("No filename was given to save the image to.");
        if (error_is_fatal) {
            exit(EXIT_FAILURE);
        }
        return -1;
    }
    return 0;
}
