#include "exec.h"

void exec() {
    int i;
    screen _screen;
    struct stack *s;
    struct matrix *trans_mat;
    struct matrix *pts;
    color c;
    c.red = 30;
    c.blue = 100;
    c.green = 155;
    #ifdef DEBUG
    print_debug("Num ops: %d", lastop);
    #endif
    _screen = new_screen();
    s = new_stack();
    trans_mat = NULL;
    pts = new_matrix(4, 4);
    for (i = 0; i < lastop; ++i) {
        struct command current_op = op[i];
        struct matrix *tmp;
        switch (current_op.opcode) {
            case PUSH:
                #ifdef DEBUG
                print_debug("Got push command");
                #endif
                // Add a new matrix to the origin stack
                push(s);
                break;
            case POP:
                #ifdef DEBUG
                print_debug("Got pop command");
                #endif
                // Remove the top matrix in the origin stack
                pop(s);
                break;
            case MOVE:
                #ifdef DEBUG
                print_debug("Got move command");
                #endif
                free_matrix(trans_mat);
                // Create a transformation matrix using the data stored for the
                // translate command
                trans_mat = make_translate(current_op.op.move.d[0],
                                           current_op.op.move.d[1],
                                           current_op.op.move.d[2]);
                // Apply the transformation to the top matrix in the origin
                // stack
                tmp = matrix_mult(s->data[s->top], trans_mat);
                free_matrix(s->data[s->top]);
                free_matrix(trans_mat);
                s->data[s->top] = tmp;
                trans_mat = NULL;
                break;
            case ROTATE:
                #ifdef DEBUG
                print_debug("Got rotate command");
                #endif
                free_matrix(trans_mat);
                switch (current_op.op.rotate.axis) {
                    case X_AXIS:
                        trans_mat = make_rotX(current_op.op.rotate.degrees);
                        tmp = matrix_mult(s->data[s->top], trans_mat);
                        free_matrix(s->data[s->top]);
                        free_matrix(trans_mat);
                        s->data[s->top] = tmp;
                        break;
                    case Y_AXIS:
                        trans_mat = make_rotY(current_op.op.rotate.degrees);
                        tmp = matrix_mult(s->data[s->top], trans_mat);
                        free_matrix(s->data[s->top]);
                        free_matrix(trans_mat);
                        s->data[s->top] = tmp;
                        break;

                    case Z_AXIS:
                        trans_mat = make_rotZ(current_op.op.rotate.degrees);
                        tmp = matrix_mult(s->data[s->top], trans_mat);
                        free_matrix(s->data[s->top]);
                        free_matrix(trans_mat);
                        s->data[s->top] = tmp;
                        break;
                }
                trans_mat = NULL;
                break;
            case SCALE:
                #ifdef DEBUG
                print_debug("Got scale command");
                #endif
                free_matrix(trans_mat);
                trans_mat = make_scale(current_op.op.scale.d[0],
                                       current_op.op.scale.d[1],
                                       current_op.op.scale.d[2]);
                tmp = matrix_mult(s->data[s->top], trans_mat);
                free_matrix(s->data[s->top]);
                free_matrix(trans_mat);
                s->data[s->top] = tmp;
                trans_mat = NULL;
                break;
            case BOX:
                #ifdef DEBUG
                print_debug("Got box command");
                #endif
                add_box(pts,
                        current_op.op.box.d0[0],
                        current_op.op.box.d0[1],
                        current_op.op.box.d0[2],
                        current_op.op.box.d1[0],
                        current_op.op.box.d1[1],
                        current_op.op.box.d1[2],
                        global_draw_mode);
                tmp = matrix_mult(s->data[s->top], pts);
                draw(_screen, tmp, c);
                free_matrix(tmp);
                free_matrix(pts);
                pts = new_matrix(4, 4);
                break;
            case SPHERE:
                #ifdef DEBUG
                print_debug("Got sphere command");
                #endif
                add_sphere(pts,
                           current_op.op.sphere.step_size,
                           current_op.op.sphere.d[0],
                           current_op.op.sphere.d[1],
                           current_op.op.sphere.d[2],
                           current_op.op.sphere.r,
                           global_draw_mode);
                tmp = matrix_mult(s->data[s->top], pts);
                draw(_screen, tmp, c);
                free_matrix(tmp);
                free_matrix(pts);
                pts = new_matrix(4, 4);
                break;
            case TORUS:
                #ifdef DEBUG
                print_debug("Got torus command");
                #endif
                add_torus(pts,
                          current_op.op.torus.step_size,
                          current_op.op.torus.d[0],
                          current_op.op.torus.d[1],
                          current_op.op.torus.d[2],
                          current_op.op.torus.circle_radius,
                          current_op.op.torus.torus_radius,
                          global_draw_mode);
                tmp = matrix_mult(s->data[s->top], pts);
                draw(_screen, tmp, c);
                free_matrix(tmp);
                free_matrix(pts);
                pts = new_matrix(4, 4);
                break;
            case LINE:
                #ifdef DEBUG
                print_debug("Got line command");
                #endif
                add_edge(pts,
                         current_op.op.line.p0[0],
                         current_op.op.line.p0[1],
                         current_op.op.line.p0[2],
                         current_op.op.line.p1[0],
                         current_op.op.line.p1[1],
                         current_op.op.line.p1[2]);
                tmp = matrix_mult(s->data[s->top], pts);
                draw(_screen, tmp, c);
                free_matrix(tmp);
                free_matrix(pts);
                pts = new_matrix(4, 4);
                break;
            case SAVE:
                #ifdef DEBUG
                print_debug("Got save command");
                #endif
                save_extension(_screen, current_op.op.save.p->name);
                break;
            case DISPLAY:
                #ifdef DEBUG
                print_debug("Got display command");
                #endif
                display(_screen);
                break;
            case DRAW_MODE: ; // Obligatory empty statement
                #ifdef DEBUG
                print_debug("Got draw-mode command");
                #endif
                char *mode = current_op.op.drawmode.p->name;
                if (strcmp(mode, "lines") == 0) {
                    global_draw_mode = DRAW_LINE;
                }
                else if (strcmp(mode, "polygons") == 0) {
                    global_draw_mode = DRAW_POLYGON;
                }
                else {
                    print_warning("Invalid argument for plot mode command: \"%s\"", mode);
                }
                break;
            case RESIZE:
                #ifdef DEBUG
                print_debug("Got resize command");
                #endif
                _screen = resize_screen(_screen,
                                        current_op.op.resize.x,
                                        current_op.op.resize.y);
                break;

            default:
                break;
        }
    }
    free_screen(_screen);
    free_stack(s);
    free_matrix(trans_mat);
    free_matrix(pts);
}

static void draw(screen s, struct matrix *pts, color c) {
    // Draw the points matrix to the screen using the current drawing mode
    if (global_plot_mode == PLOT_CARTESIAN) {
        color axis_color;
        axis_color.red = 255;
        axis_color.blue = 0;
        axis_color.green = 0;
        draw_axes(s, axis_color);
    }
    switch (global_draw_mode) {
        case DRAW_LINE:
            draw_lines(s, c, pts, global_plot_mode);
            break;
        case DRAW_POLYGON:
            draw_polygons(s, c, pts, global_plot_mode);
            break;
    }
}
