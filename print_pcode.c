#include "print_pcode.h"

int lastop;
struct command op[MAX_COMMANDS];

void print_pcode() {
    int i;
    for (i=0;i<lastop;i++) {
        printf("%d: ",i);
        struct command current_op = op[i];
        switch (current_op.opcode) {
            case LIGHT:
                printf("Light: %s at: %6.2f %6.2f %6.2f",
                        current_op.op.light.p->name,
                        current_op.op.light.c[0], current_op.op.light.c[1],
                        current_op.op.light.c[2]);
                break;
            case AMBIENT:
                printf("Ambient: %6.2f %6.2f %6.2f",
                        current_op.op.ambient.c[0],
                        current_op.op.ambient.c[1],
                        current_op.op.ambient.c[2]);
                break;

            case CONSTANTS:
                printf("Constants: %s", current_op.op.constants.p->name);
                break;
            case SAVE_COORDS:
                printf("Save Coords: %s", current_op.op.save_coordinate_system.p->name);
                break;
            case CAMERA:
                printf("Camera: eye: %6.2f %6.2f %6.2f\taim: %6.2f %6.2f %6.2f",
                        current_op.op.camera.eye[0], current_op.op.camera.eye[1],
                        current_op.op.camera.eye[2],
                        current_op.op.camera.aim[0], current_op.op.camera.aim[1],
                        current_op.op.camera.aim[2]);

                break;
            case SPHERE:
                printf("Sphere: x=%6.2f y=%6.2f z=%6.2f r=%6.2f step=%6.2f",
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
                printf("Torus: x=%6.2f y=%6.2f z=%6.2f"
                       "circle_radius=%6.2f torus_radius=%6.2f step=%6.2f",
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
                printf("Box: d0: x0=%6.2f y0=%6.2f z0=%6.2f d1: x1=%6.2f y1=%6.2f z1=%6.2f",
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
                printf("Line: from: %6.2f %6.2f %6.2f to: %6.2f %6.2f %6.2f",
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
                printf("Mesh: filename: %s", current_op.op.mesh.name);
                if (current_op.op.mesh.constants != NULL) {
                    printf("\tconstants: %s", current_op.op.mesh.constants->name);
                }
                break;
            case SET:
                printf("Set: %s %6.2f",
                        current_op.op.set.p->name,
                        current_op.op.set.p->s.value);
                break;
            case MOVE:
                printf("Move: x=%6.2f y=%6.2f z=%6.2f",
                        current_op.op.move.d[0], current_op.op.move.d[1],
                        current_op.op.move.d[2]);
                if (current_op.op.move.p != NULL) {
                    printf("\tknob: %s", current_op.op.move.p->name);
                }
                break;
            case SCALE:
                printf("Scale: %6.2f %6.2f %6.2f",
                        current_op.op.scale.d[0], current_op.op.scale.d[1],
                        current_op.op.scale.d[2]);
                if (current_op.op.scale.p != NULL) {
                    printf("\tknob: %s", current_op.op.scale.p->name);
                }
                break;
            case ROTATE:
                printf("Rotate: axis: %d degrees: %6.2f",
                        current_op.op.rotate.axis,
                        current_op.op.rotate.degrees);
                if (current_op.op.rotate.p != NULL) {
                    printf("\tknob: %s", current_op.op.rotate.p->name);
                }
                break;
            case BASENAME:
                printf("Basename: %s", current_op.op.basename.p->name);
                break;
            case SAVE_KNOBS:
                printf("Save knobs: %s", current_op.op.save_knobs.p->name);
                break;
            case TWEEN:
                printf("Tween: %4.0f %4.0f, %s %s",
                        current_op.op.tween.start_frame,
                        current_op.op.tween.end_frame,
                        current_op.op.tween.knob_list0->name,
                        current_op.op.tween.knob_list1->name);
                break;
            case FRAMES:
                printf("Num frames: %4.0f", current_op.op.frames.num_frames);
                break;
            case VARY:
                printf("Vary: %4.0f %4.0f, %4.0f %4.0f",
                        current_op.op.vary.start_frame,
                        current_op.op.vary.end_frame,
                        current_op.op.vary.start_val,
                        current_op.op.vary.end_val);
                break;
            case PUSH:
                printf("Push");
                break;
            case POP:
                printf("Pop");
                break;
            case GENERATE_RAYFILES:
                printf("Generate Ray Files");
                break;
            case SAVE:
                printf("Save: %s", current_op.op.save.p->name);
                break;
            case SHADING:
                printf("Shading: %s", current_op.op.shading.p->name);
                break;
            case SETKNOBS:
                printf("Setknobs: %f", current_op.op.setknobs.value);
                break;
            case FOCAL:
                printf("Focal: %f", current_op.op.focal.value);
                break;
            case DISPLAY:
                printf("Display");
                break;
            case DRAW_MODE:
                printf("Draw mode: %s", current_op.op.drawmode.p->name);
                break;
            case RESIZE:
                printf("Resize: x=%d y=%d", current_op.op.resize.x,
                                            current_op.op.resize.y);
        }
        printf("\n");
    }
}


