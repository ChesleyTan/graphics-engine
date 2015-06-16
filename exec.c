#include "exec.h"

screen _screen;
int num_frames = 0;
char *default_basename = "frame";
char *basename = NULL;
struct vary_node **vary_knobs;

screen exec(char return_screen) {
    int i;
    struct stack *s;
    struct matrix *trans_mat;
    struct matrix *pts;
    color c;
    allocate_z_buffer();
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
                if (current_op.op.move.p != NULL) { // If a knob exists, use it
                    // Only apply translation to coordinates axes who initial values
                    // are non-zero
                    trans_mat = make_translate(current_op.op.move.d[0] == 0 ? 0
                                               : current_op.op.move.p->s.value,
                                               current_op.op.move.d[1] == 0 ? 0
                                               : current_op.op.move.p->s.value,
                                               current_op.op.move.d[2] == 0 ? 0
                                               : current_op.op.move.p->s.value);
                }
                else {
                    trans_mat = make_translate(current_op.op.move.d[0],
                                               current_op.op.move.d[1],
                                               current_op.op.move.d[2]);
                }
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
                        if (current_op.op.rotate.p != NULL) { // If knob exists, use knob value
                            trans_mat = make_rotX(current_op.op.rotate.p->s.value);
                        }
                        else {
                            trans_mat = make_rotX(current_op.op.rotate.degrees);
                        }
                        break;
                    case Y_AXIS:
                        if (current_op.op.rotate.p != NULL) { // If knob exists, use knob value
                            trans_mat = make_rotY(current_op.op.rotate.p->s.value);
                        }
                        else {
                            trans_mat = make_rotY(current_op.op.rotate.degrees);
                        }
                        break;

                    case Z_AXIS:
                        if (current_op.op.rotate.p != NULL) { // If knob exists, use knob value
                            trans_mat = make_rotZ(current_op.op.rotate.p->s.value);
                        }
                        else {
                            trans_mat = make_rotZ(current_op.op.rotate.degrees);
                        }
                        break;
                }
                tmp = matrix_mult(s->data[s->top], trans_mat);
                free_matrix(s->data[s->top]);
                free_matrix(trans_mat);
                s->data[s->top] = tmp;
                trans_mat = NULL;
                break;
            case SCALE:
                #ifdef DEBUG
                print_debug("Got scale command");
                #endif
                free_matrix(trans_mat);
                if (current_op.op.scale.p != NULL) { // If knob exists, use it
                    // Only apply scale to coordinates axes who initial values
                    // are non-zero
                    trans_mat = make_scale(current_op.op.scale.d[0] == 0 ? 1
                                           : current_op.op.scale.p->s.value,
                                           current_op.op.scale.d[1] == 0 ? 1
                                           : current_op.op.scale.p->s.value,
                                           current_op.op.scale.d[2] == 0 ? 1
                                           : current_op.op.scale.p->s.value);
                }
                else {
                    trans_mat = make_scale(current_op.op.scale.d[0],
                                           current_op.op.scale.d[1],
                                           current_op.op.scale.d[2]);
                }
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
                char *draw_mode = current_op.op.drawmode.p->name;
                if (strcmp(draw_mode, "lines") == 0) {
                    global_draw_mode = DRAW_LINE;
                }
                else if (strcmp(draw_mode, "polygons") == 0) {
                    global_draw_mode = DRAW_POLYGON;
                }
                else {
                    print_error("Invalid argument for draw-mode command: \"%s\"", draw_mode);
                }
                break;
            case RENDER_MODE: ; // Obligatory empty statement
                #ifdef DEBUG
                print_debug("Got render-mode command");
                #endif
                char *render_mode = current_op.op.rendermode.p->name;
                if (strcmp(render_mode, "wireframe") == 0) {
                    global_render_mode = RENDER_WIREFRAME;
                }
                else if (strcmp(render_mode, "surface") == 0) {
                    global_render_mode = RENDER_SURFACE;
                }
                else {
                    print_error("Invalid argument for render-mode command: \"%s\"", render_mode);
                }
                break;
            case SHADE_MODE: ; // Obligatory empty statement
                #ifdef DEBUG
                print_debug("Got shade-mode command");
                #endif
                char *shade_mode = current_op.op.drawmode.p->name;
                if (strcmp(shade_mode, "flat") == 0) {
                    global_shade_mode = SHADE_FLAT;
                }
                else if (strcmp(shade_mode, "gouraud") == 0) {
                    global_shade_mode = SHADE_GOURAUD;
                }
                else if (strcmp(shade_mode, "phong") == 0) {
                    global_shade_mode = SHADE_PHONG;
                }
                else {
                    print_error("Invalid argument for shade-mode command: \"%s\"", shade_mode);
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
            case LIGHT:
                #ifdef DEBUG
                print_debug("Got light command");
                #endif
                set_light_vector(current_op.op.light.coord[0],
                                current_op.op.light.coord[1],
                                current_op.op.light.coord[2]);
                break;
            case CONSTANTS:
                #ifdef DEBUG
                print_debug("Got constants command");
                #endif
                set_lighting_constants(current_op.op.constants.p->s.c);
                break;
            case FRAMES: // handled by parse_animation_cmds()
                break;
            case BASENAME: // handled by parse_animation_cmds()
                break;
            case VARY: // handled by parse_animation_cmds()
                break;
            default:
                print_error("This command has not been implemented yet.");
                free_stack(s);
                free_matrix(trans_mat);
                free_matrix(pts);
                free_all();
                exit(EXIT_FAILURE);
                break;
        }
    }
    free_stack(s);
    free_matrix(trans_mat);
    free_matrix(pts);
    free_z_buffer();
    if (!return_screen) {
        free_screen(_screen);
        return NULL;
    }
    else {
        return _screen;
    }
}

void parse_animation_cmds() {
    int i, u;
    for (i = 0; i < lastop; ++i) {
        struct command current_op = op[i];
        switch (current_op.opcode) {
            case FRAMES:
                num_frames = current_op.op.frames.num_frames;
                #ifdef DEBUG
                print_debug("num_frames: %d", num_frames);
                #endif
                allocate_vary_knobs();
                break;
            case BASENAME:
                basename = current_op.op.basename.p->name;
                #ifdef DEBUG
                print_debug("basename: %s", basename);
                #endif
                break;
            case VARY: ; // Obligatory empty statement
                double vary_slope = (current_op.op.vary.end_val
                                  - current_op.op.vary.start_val)
                                    / (current_op.op.vary.end_frame
                                       - current_op.op.vary.start_frame);
                // Populate the vary_knobs array
                for (u = 0; u < num_frames; ++u) {
                    if (u < current_op.op.vary.start_frame) {
                        // Don't overwrite previously set values
                        if (vary_node_uniq(u, current_op.op.vary.p)) {
                            // Add a vary_node with value of the start_val
                            struct vary_node *new = (struct vary_node *)
                                       malloc(sizeof(struct vary_node));
                            new->knob = current_op.op.vary.p;
                            new->next = NULL;
                            new->next_value = current_op.op.vary.start_val;
                            get_vary_knobs_tail(u)->next = new;
                        }
                    }
                    else if (u > current_op.op.vary.end_frame) {
                        // Don't overwrite previously set values
                        if (vary_node_uniq(u, current_op.op.vary.p)) {
                            // Add a vary_node with value of the end_val
                            struct vary_node *new = (struct vary_node *)
                                       malloc(sizeof(struct vary_node));
                            new->knob = current_op.op.vary.p;
                            new->next = NULL;
                            new->next_value = current_op.op.vary.end_val;
                            get_vary_knobs_tail(u)->next = new;
                        }
                    }
                    else {
                        // Perform calculation for the value
                        struct vary_node *new = (struct vary_node *)
                                    malloc(sizeof(struct vary_node));
                        new->knob = current_op.op.vary.p;
                        new->next = NULL;
                        new->next_value = current_op.op.vary.start_val
                                        + (u - current_op.op.vary.start_frame)
                                          * vary_slope;
                        get_vary_knobs_tail(u)->next = new;
                    }
                }
                break;
        }
    }
    // Make sure a basename is set; If a basename is not set, then use a
    // default value.
    if (basename == NULL) {
        basename = default_basename;
        print_warning("The basename for generated frame files was not specified! "
                      "The default value \"%s\" will be used.", default_basename);
    }
    #ifdef DEBUG
    for (i = 0; i < num_frames; ++i) {
        print_debug("Frame: %d", i);
        struct vary_node *curr = vary_knobs[i]->next;
        while (curr) {
            print_debug("\tKnob: %s; val: %lf",
                   curr->knob->name,
                   curr->next_value);
            curr = curr->next;
        }
    }
    #endif
}

void exec_animation() {
    screen s;
    int frame;
    char save_filename[256];
    if (access("frames/", F_OK) != 0) {
        // Directory does not exist, so create it
        if (errno == ENOENT) {
            int ret = mkdir("frames", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); // 775
            if (ret != 0) {
                print_errno("Error encountered while creating `frames/` directory.");
                free_all();
                exit(EXIT_FAILURE);
            }
        }
        // If file is not a directory, print error
        else if (errno == ENOENT) {
            print_error("`frames` already exists, but is not a directory! Cannot save animation output. Exiting....");
            free_all();
            exit(EXIT_FAILURE);
        }
    }
    for (frame = 0; frame < num_frames; ++frame) {
        // Set the knobs to their next value for the current frame
        struct vary_node *n = vary_knobs[frame]->next;
        while (n != NULL) {
            n->knob->s.value = n->next_value;
            #ifdef DEBUG
            print_debug("Knob: %s", n->knob->name);
            print_debug("Value: %lf", n->knob->s.value);
            #endif
            n = n->next;
        }
        // Execute opcodes
        s = exec(TRUE);
        // Save screen to file
        snprintf(save_filename, sizeof(save_filename), "frames/%s%03d.png", basename, frame);
        save_filename[sizeof(save_filename) - 1] = '\0';
        save_extension(s, save_filename);
        free_screen(s);
    }
    free_vary_knobs();
    return;
}

struct vary_node *get_vary_knobs_tail(int frame) {
    struct vary_node *curr = vary_knobs[frame];
    while (curr->next) {
        curr = curr->next;
    }
    return curr;
}

char vary_node_uniq(int frame, SYMTAB *knob) {
    struct vary_node *curr = vary_knobs[frame];
    while (curr->next) {
        curr = curr->next;
        if (curr->knob == knob) {
            return FALSE;
        }
    }
    return TRUE;
}

void allocate_vary_knobs() {
    // Allocate memory for vary_knobs array
    int i;
    vary_knobs = (struct vary_node **) malloc(
                    num_frames * sizeof(struct vary_node *));
    // Initialize all pointers of the array with a dummy node
    for (i = 0; i < num_frames; ++i) {
        vary_knobs[i] = (struct vary_node *) malloc(sizeof(struct vary_node));
        vary_knobs[i]->knob = NULL;
        vary_knobs[i]->next_value = 0;
        vary_knobs[i]->next = NULL;
    }
}

void free_exec_screen() {
    if (_screen != NULL) {
        free_screen(_screen);
        _screen = NULL;
    }
}

void free_vary_knobs() {
    if (vary_knobs != NULL) {
        #ifdef DEBUG
        print_debug("Freeing vary_knobs");
        #endif
        int i;
        for (i = 0; i < num_frames; ++i) {
            struct vary_node *curr = vary_knobs[i];
            struct vary_node *next = curr->next;
            while (next) {
                free(curr);
                curr = next;
                next = curr->next;
            }
            free(curr);
        }
        free(vary_knobs);
        vary_knobs = NULL;
    }
    return;
}

void free_all() {
    // TODO update this when new global dynamically allocated blocks are
    // added
    #ifdef DEBUG
    print_debug("Calling free_all()");
    #endif
    free_vary_knobs();
    free_exec_screen();
    free_z_buffer();
    free_table();
}
