#include "print_pcode.h"

int lastop;
struct command op[MAX_COMMANDS];

int display_length(int n) {
    int l = 0;
    while (n > 0) {
        n /= 10;
        ++l;
    }
    return l;
}

void print_pcode() {
    int i;
    int max_lineno_length = display_length(lastop);
    for (i=0;i<lastop;i++) {
        printf("%s%s%s%s %-*d %s ", bold_prefix, fg_green_118, prefix, bg_black, max_lineno_length, i+1, reset);
        struct command current_op = op[i];
        switch (current_op.opcode) {
            case LIGHT:
                printf("%s%sLight%s: at %6.2f %6.2f %6.2f",
                        bold_prefix, fg_blue_30, reset,
                        current_op.op.light.coord[0], current_op.op.light.coord[1],
                        current_op.op.light.coord[2]);
                break;
            case AMBIENT:
                printf("%s%sAmbient%s: %6.2f %6.2f %6.2f",
                        bold_prefix, fg_blue_30, reset,
                        current_op.op.ambient.c[0],
                        current_op.op.ambient.c[1],
                        current_op.op.ambient.c[2]);
                break;

            case CONSTANTS:
                printf("%s%sConstants%s: %s\n",
                        bold_prefix, fg_blue_30, reset,
                        current_op.op.constants.p->name);
                print_constants(current_op.op.constants.p->s.c);
                break;
            case SAVE_COORDS:
                printf("%s%sSave Coords%s: %s",
                        bold_prefix, fg_blue_30, reset,
                        current_op.op.save_coordinate_system.p->name);
                break;
            case CAMERA:
                printf("%s%sCamera%s: eye: %6.2f %6.2f %6.2f\taim: %6.2f %6.2f %6.2f",
                        bold_prefix, fg_blue_30, reset,
                        current_op.op.camera.eye[0], current_op.op.camera.eye[1],
                        current_op.op.camera.eye[2],
                        current_op.op.camera.aim[0], current_op.op.camera.aim[1],
                        current_op.op.camera.aim[2]);

                break;
            case SPHERE:
                printf("%s%sSphere%s: x=%6.2f y=%6.2f z=%6.2f r=%6.2f step=%6.2f",
                        bold_prefix, fg_blue_30, reset,
                        current_op.op.sphere.d[0], current_op.op.sphere.d[1],
                        current_op.op.sphere.d[2],
                        current_op.op.sphere.r,
                        current_op.op.sphere.step_size);
                if (current_op.op.sphere.constants != NULL) {
                    printf("\tconstants: %s", current_op.op.sphere.constants->name);
                }
                if (current_op.op.sphere.cs != NULL) {
                    printf("\tcs: %s", current_op.op.sphere.cs->name);
                }

                break;
            case TORUS:
                printf("%s%sTorus%s: x=%6.2f y=%6.2f z=%6.2f"
                       "circle_radius=%6.2f torus_radius=%6.2f step=%6.2f",
                        bold_prefix, fg_blue_30, reset,
                        current_op.op.torus.d[0], current_op.op.torus.d[1],
                        current_op.op.torus.d[2],
                        current_op.op.torus.circle_radius,
                        current_op.op.torus.torus_radius,
                        current_op.op.torus.step_size);
                if (current_op.op.torus.constants != NULL)
                {
                    printf("\tconstants: %s", current_op.op.torus.constants->name);
                }
                if (current_op.op.torus.cs != NULL)
                {
                    printf("\tcs: %s", current_op.op.torus.cs->name);
                }

                break;
            case BOX:
                printf("%s%sBox%s: d0: x0=%6.2f y0=%6.2f z0=%6.2f d1: x1=%6.2f y1=%6.2f z1=%6.2f",
                        bold_prefix, fg_blue_30, reset,
                        current_op.op.box.d0[0], current_op.op.box.d0[1],
                        current_op.op.box.d0[2],
                        current_op.op.box.d1[0], current_op.op.box.d1[1],
                        current_op.op.box.d1[2]);
                if (current_op.op.box.constants != NULL)
                {
                    printf("\tconstants: %s", current_op.op.box.constants->name);
                }
                if (current_op.op.box.cs != NULL)
                {
                    printf("\tcs: %s", current_op.op.box.cs->name);
                }

                break;
            case LINE:
                printf("%s%sLine%s: from: %6.2f %6.2f %6.2f to: %6.2f %6.2f %6.2f",
                        bold_prefix, fg_blue_30, reset,
                        current_op.op.line.p0[0], current_op.op.line.p0[1],
                        current_op.op.line.p0[1],
                        current_op.op.line.p1[0], current_op.op.line.p1[1],
                        current_op.op.line.p1[1]);
                if (current_op.op.line.constants != NULL) {
                    printf("\n\tConstants: %s", current_op.op.line.constants->name);
                }
                if (current_op.op.line.cs0 != NULL) {
                    printf("\n\tCS0: %s", current_op.op.line.cs0->name);
                }
                if (current_op.op.line.cs1 != NULL) {
                    printf("\n\tCS1: %s", current_op.op.line.cs1->name);
                }
                break;
            case MESH:
                printf("%s%sMesh%s: filename: %s",
                        bold_prefix, fg_blue_30, reset,
                        current_op.op.mesh.name);
                if (current_op.op.mesh.constants != NULL) {
                    printf("\tconstants: %s", current_op.op.mesh.constants->name);
                }
                break;
            case SET:
                printf("%s%sSet%s: %s %6.2f",
                        bold_prefix, fg_blue_30, reset,
                        current_op.op.set.p->name,
                        current_op.op.set.p->s.value);
                break;
            case MOVE:
                printf("%s%sMove%s: x=%6.2f y=%6.2f z=%6.2f",
                        bold_prefix, fg_blue_30, reset,
                        current_op.op.move.d[0], current_op.op.move.d[1],
                        current_op.op.move.d[2]);
                if (current_op.op.move.p != NULL) {
                    printf("\tknob: %s", current_op.op.move.p->name);
                }
                break;
            case SCALE:
                printf("%s%sScale%s: %6.2f %6.2f %6.2f",
                        bold_prefix, fg_blue_30, reset,
                        current_op.op.scale.d[0], current_op.op.scale.d[1],
                        current_op.op.scale.d[2]);
                if (current_op.op.scale.p != NULL) {
                    printf("\tknob: %s", current_op.op.scale.p->name);
                }
                break;
            case ROTATE:
                printf("%s%sRotate%s: axis: %d degrees: %6.2f",
                        bold_prefix, fg_blue_30, reset,
                        current_op.op.rotate.axis,
                        current_op.op.rotate.degrees);
                if (current_op.op.rotate.p != NULL) {
                    printf("\tknob: %s", current_op.op.rotate.p->name);
                }
                break;
            case BASENAME:
                printf("%s%sBasename%s: %s",
                        bold_prefix, fg_blue_30, reset,
                        current_op.op.basename.p->name);
                break;
            case SAVE_KNOBS:
                printf("%s%sSave knobs%s: %s",
                        bold_prefix, fg_blue_30, reset,
                        current_op.op.save_knobs.p->name);
                break;
            case TWEEN:
                printf("%s%sTween%s: %4.0f %4.0f, %s %s",
                        bold_prefix, fg_blue_30, reset,
                        current_op.op.tween.start_frame,
                        current_op.op.tween.end_frame,
                        current_op.op.tween.knob_list0->name,
                        current_op.op.tween.knob_list1->name);
                break;
            case FRAMES:
                printf("%s%sNum frames%s: %d",
                        bold_prefix, fg_blue_30, reset,
                        current_op.op.frames.num_frames);
                break;
            case VARY:
                printf("%s%sVary%s: %d %d, %4.0f %4.0f",
                        bold_prefix, fg_blue_30, reset,
                        current_op.op.vary.start_frame,
                        current_op.op.vary.end_frame,
                        current_op.op.vary.start_val,
                        current_op.op.vary.end_val);
                break;
            case PUSH:
                printf("%s%sPush%s", bold_prefix, fg_blue_30, reset);
                break;
            case POP:
                printf("%s%sPop%s", bold_prefix, fg_blue_30, reset);
                break;
            case GENERATE_RAYFILES:
                printf("%s%sGenerate Ray Files%s", bold_prefix, fg_blue_30, reset);
                break;
            case SAVE:
                printf("%s%sSave%s: %s",
                        bold_prefix, fg_blue_30, reset,
                        current_op.op.save.p->name);
                break;
            case SETKNOBS:
                printf("%s%sSetknobs%s: %f",
                        bold_prefix, fg_blue_30, reset,
                        current_op.op.setknobs.value);
                break;
            case FOCAL:
                printf("%s%sFocal%s: %f",
                        bold_prefix, fg_blue_30, reset,
                        current_op.op.focal.value);
                break;
            case DISPLAY:
                printf("%s%sDisplay%s", bold_prefix, fg_blue_30, reset);
                break;
            case DRAW_MODE:
                printf("%s%sDraw mode%s: %s",
                        bold_prefix, fg_blue_30, reset,
                        current_op.op.drawmode.p->name);
                break;
            case RENDER_MODE:
                printf("%s%sRender mode%s: %s",
                        bold_prefix, fg_blue_30, reset,
                        current_op.op.rendermode.p->name);
                break;
            case SHADE_MODE:
                printf("%s%sShade mode%s: %s",
                        bold_prefix, fg_blue_30, reset,
                        current_op.op.shading.p->name);
                break;
            case RESIZE:
                printf("%s%sResize%s: x=%d y=%d",
                        bold_prefix, fg_blue_30, reset,
                        current_op.op.resize.x,
                        current_op.op.resize.y);
        }
        printf("\n");
    }
}


