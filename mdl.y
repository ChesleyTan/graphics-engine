%{
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include "parser.h"
#include "matrix.h"
#include "exec.h"
#include "utils.h"

#define YYERROR_VERBOSE 1
#ifdef DEBUG
#define YYDEBUG 1
#endif

SYMTAB *s;
struct light *l;
struct constants *c;
struct command op[MAX_COMMANDS];
struct matrix *m;
int lastop = 0;
int lineno = 0;
char is_animation = FALSE;

// Bison headers and variables
void yyerror(const char *);
extern FILE *yyin;

%}

/* =============== DEFINITIONS ============= */

%union {
    double val;
    char string[255];
}

%token COMMENT
%token <val> NUMBER
%token <string> DRAW_MODE DRAW_TYPE RENDER_MODE RENDER_TYPE RESIZE
%token <string> LIGHT AMBIENT
%token <string> CONSTANTS SAVE_COORDS CAMERA 
%token <string> SPHERE TORUS BOX LINE CS MESH TEXTURE
%token <string> STRING
%token <string> SET MOVE SCALE ROTATE BASENAME SAVE_KNOBS TWEEN FRAMES VARY 
%token <string> PUSH POP SAVE GENERATE_RAYFILES
%token <string> SHADE_MODE SHADING_TYPE SETKNOBS FOCAL DISPLAY WEB
%token <string> CO
%%
/* ============= END DEFINITIONS ============= */

/* ================== RULES ================== */

input:
| input command
;

command: 
  COMMENT {
  ++lineno;
  }

| DRAW_MODE DRAW_TYPE {
  ++lineno;
  op[lastop].opcode=DRAW_MODE;
  op[lastop].op.drawmode.p = add_symbol($2,SYM_STRING,NULL);
  ++lastop;
}

| RENDER_MODE RENDER_TYPE {
  ++lineno;
  op[lastop].opcode=RENDER_MODE;
  op[lastop].op.rendermode.p = add_symbol($2,SYM_STRING,NULL);
  ++lastop;
}

| SHADE_MODE SHADING_TYPE {
  ++lineno;
  op[lastop].opcode = SHADE_MODE;
  op[lastop].op.shading.p = add_symbol($2,SYM_STRING,NULL);
  ++lastop;
}

| RESIZE NUMBER NUMBER {
  ++lineno;
  op[lastop].opcode=RESIZE;
  op[lastop].op.resize.x = round($2);
  op[lastop].op.resize.y = round($3);
  ++lastop;
}

    /* TODO allow knobs to set light position */
| LIGHT NUMBER NUMBER NUMBER {
  ++lineno;
  op[lastop].opcode=LIGHT;
  op[lastop].op.light.coord[0] = $2;
  op[lastop].op.light.coord[1] = $3;
  op[lastop].op.light.coord[2] = $4;
  ++lastop;
}

| MOVE NUMBER NUMBER NUMBER STRING { 
  ++lineno;
  op[lastop].opcode = MOVE;
  op[lastop].op.move.d[0] = $2;
  op[lastop].op.move.d[1] = $3;
  op[lastop].op.move.d[2] = $4;
  op[lastop].op.move.d[3] = 0;
  op[lastop].op.move.p = add_symbol($5,SYM_STRING,NULL);
  ++lastop;
}

| MOVE NUMBER NUMBER NUMBER {
  ++lineno;
  op[lastop].opcode = MOVE;
  op[lastop].op.move.d[0] = $2;
  op[lastop].op.move.d[1] = $3;
  op[lastop].op.move.d[2] = $4;
  op[lastop].op.move.d[3] = 0;
  op[lastop].op.move.p = NULL;
  ++lastop;
}

| CONSTANTS STRING NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER {
  ++lineno;
  c = (struct constants *)malloc(sizeof(struct constants));
  c->iar = (int) $3;
  c->iag = (int) $4;
  c->iab = (int) $5;
  c->kar = $6;
  c->kag = $7;
  c->kab = $8;
  c->idr = (int) $9;
  c->idg = (int) $10;
  c->idb = (int) $11;
  c->kdr = $12;
  c->kdg = $13;
  c->kdb = $14;
  c->isr = (int) $15;
  c->isg = (int) $16;
  c->isb = (int) $17;
  c->ksr = $18;
  c->ksg = $19;
  c->ksb = $20;
  c->spec_expt = (int) $21;
  op[lastop].op.constants.p =  add_symbol($2,SYM_CONSTANTS,c);
  op[lastop].opcode=CONSTANTS;
  ++lastop;
}

| SAVE_COORDS STRING {
  ++lineno;
  op[lastop].opcode = SAVE_COORDS;
  m = new_matrix(4,4);
  op[lastop].op.save_coordinate_system.p = add_symbol($2,SYM_MATRIX,m);
  ++lastop;
}

| CAMERA NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER {
  ++lineno;
  op[lastop].opcode = CAMERA;
  op[lastop].op.camera.eye[0] = $2;
  op[lastop].op.camera.eye[1] = $3;
  op[lastop].op.camera.eye[2] = $4;
  op[lastop].op.camera.eye[3] = 0;
  op[lastop].op.camera.aim[0] = $5;
  op[lastop].op.camera.aim[1] = $6;
  op[lastop].op.camera.aim[2] = $7;
  op[lastop].op.camera.aim[3] = 0;
  ++lastop;
}

| TEXTURE STRING NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER {
  ++lineno;
  op[lastop].opcode = TEXTURE;
  c = (struct constants *)malloc(sizeof(struct constants));
  op[lastop].op.texture.d0[0] = $3;
  op[lastop].op.texture.d0[1] = $4;
  op[lastop].op.texture.d0[2] = $5;
  op[lastop].op.texture.d1[0] = $6;
  op[lastop].op.texture.d1[1] = $7;
  op[lastop].op.texture.d1[2] = $8;
  op[lastop].op.texture.d2[0] = $9;
  op[lastop].op.texture.d2[1] = $10;
  op[lastop].op.texture.d2[2] = $11;
  op[lastop].op.texture.d3[0] = $12;
  op[lastop].op.texture.d3[1] = $13;
  op[lastop].op.texture.d3[2] = $14;
  op[lastop].op.texture.cs = NULL;
  op[lastop].op.texture.constants =  add_symbol("",SYM_CONSTANTS,c);
  op[lastop].op.texture.p = add_symbol($2,SYM_STRING,NULL);
  ++lastop;
}

| SPHERE NUMBER NUMBER NUMBER NUMBER {
  ++lineno;
  op[lastop].opcode = SPHERE;
  op[lastop].op.sphere.d[0] = $2;
  op[lastop].op.sphere.d[1] = $3;
  op[lastop].op.sphere.d[2] = $4;
  op[lastop].op.sphere.d[3] = 0;
  op[lastop].op.sphere.r = $5;
  op[lastop].op.sphere.step_size = DEFAULT_SPHERE_STEP_SIZE;
  op[lastop].op.sphere.constants = NULL;
  op[lastop].op.sphere.cs = NULL;
  ++lastop;
}

| SPHERE NUMBER NUMBER NUMBER NUMBER NUMBER {
  ++lineno;
  op[lastop].opcode = SPHERE;
  op[lastop].op.sphere.d[0] = $2;
  op[lastop].op.sphere.d[1] = $3;
  op[lastop].op.sphere.d[2] = $4;
  op[lastop].op.sphere.d[3] = 0;
  op[lastop].op.sphere.r = $5;
  op[lastop].op.sphere.step_size = $6;
  op[lastop].op.sphere.constants = NULL;
  op[lastop].op.sphere.cs = NULL;
  ++lastop;
}

| SPHERE NUMBER NUMBER NUMBER NUMBER STRING {
  ++lineno;
  op[lastop].opcode = SPHERE;
  op[lastop].op.sphere.d[0] = $2;
  op[lastop].op.sphere.d[1] = $3;
  op[lastop].op.sphere.d[2] = $4;
  op[lastop].op.sphere.d[3] = 0;
  op[lastop].op.sphere.r = $5;
  op[lastop].op.sphere.step_size = DEFAULT_SPHERE_STEP_SIZE;
  op[lastop].op.sphere.constants = NULL;
  m = new_matrix(4,4);
  op[lastop].op.sphere.cs = add_symbol($6,SYM_MATRIX,m);
  ++lastop;
}

| SPHERE NUMBER NUMBER NUMBER NUMBER NUMBER STRING {
  ++lineno;
  op[lastop].opcode = SPHERE;
  op[lastop].op.sphere.d[0] = $2;
  op[lastop].op.sphere.d[1] = $3;
  op[lastop].op.sphere.d[2] = $4;
  op[lastop].op.sphere.d[3] = 0;
  op[lastop].op.sphere.r = $5;
  op[lastop].op.sphere.step_size = $6;
  op[lastop].op.sphere.constants = NULL;
  m = new_matrix(4,4);
  op[lastop].op.sphere.cs = add_symbol($7,SYM_MATRIX,m);
  ++lastop;
}

| SPHERE STRING NUMBER NUMBER NUMBER NUMBER {
  ++lineno;
  op[lastop].opcode = SPHERE;
  op[lastop].op.sphere.d[0] = $3;
  op[lastop].op.sphere.d[1] = $4;
  op[lastop].op.sphere.d[2] = $5;
  op[lastop].op.sphere.d[3] = 0;
  op[lastop].op.sphere.r = $6;
  op[lastop].op.sphere.step_size = DEFAULT_SPHERE_STEP_SIZE;
  op[lastop].op.sphere.cs = NULL;
  c = (struct constants *)malloc(sizeof(struct constants));
  op[lastop].op.sphere.constants = add_symbol($2,SYM_CONSTANTS,c);
  ++lastop;
}

| SPHERE STRING NUMBER NUMBER NUMBER NUMBER NUMBER {
  ++lineno;
  op[lastop].opcode = SPHERE;
  op[lastop].op.sphere.d[0] = $3;
  op[lastop].op.sphere.d[1] = $4;
  op[lastop].op.sphere.d[2] = $5;
  op[lastop].op.sphere.d[3] = 0;
  op[lastop].op.sphere.r = $6;
  op[lastop].op.sphere.step_size = $7;
  op[lastop].op.sphere.cs = NULL;
  c = (struct constants *)malloc(sizeof(struct constants));
  op[lastop].op.sphere.constants = add_symbol($2,SYM_CONSTANTS,c);
  ++lastop;
}

| SPHERE STRING NUMBER NUMBER NUMBER NUMBER STRING {
  ++lineno;
  op[lastop].opcode = SPHERE;
  op[lastop].op.sphere.d[0] = $3;
  op[lastop].op.sphere.d[1] = $4;
  op[lastop].op.sphere.d[2] = $5;
  op[lastop].op.sphere.d[3] = 0;
  op[lastop].op.sphere.r = $6;
  op[lastop].op.sphere.step_size = DEFAULT_SPHERE_STEP_SIZE;
  op[lastop].op.sphere.constants = NULL;
  m = new_matrix(4,4);
  op[lastop].op.sphere.cs = add_symbol($7,SYM_MATRIX,m);
  c = (struct constants *)malloc(sizeof(struct constants));
  op[lastop].op.sphere.constants = add_symbol($2,SYM_CONSTANTS,c);
  ++lastop;
}

| SPHERE STRING NUMBER NUMBER NUMBER NUMBER NUMBER STRING {
  ++lineno;
  op[lastop].opcode = SPHERE;
  op[lastop].op.sphere.d[0] = $3;
  op[lastop].op.sphere.d[1] = $4;
  op[lastop].op.sphere.d[2] = $5;
  op[lastop].op.sphere.d[3] = 0;
  op[lastop].op.sphere.r = $6;
  op[lastop].op.sphere.step_size = $7;
  op[lastop].op.sphere.constants = NULL;
  m = new_matrix(4,4);
  op[lastop].op.sphere.cs = add_symbol($8,SYM_MATRIX,m);
  c = (struct constants *)malloc(sizeof(struct constants));
  op[lastop].op.sphere.constants = add_symbol($2,SYM_CONSTANTS,c);
  ++lastop;
}

| TORUS NUMBER NUMBER NUMBER NUMBER NUMBER {
  ++lineno;
  op[lastop].opcode = TORUS;
  op[lastop].op.torus.d[0] = $2;
  op[lastop].op.torus.d[1] = $3;
  op[lastop].op.torus.d[2] = $4;
  op[lastop].op.torus.d[3] = 0;
  op[lastop].op.torus.circle_radius = $5;
  op[lastop].op.torus.torus_radius = $6;
  op[lastop].op.torus.step_size = DEFAULT_TORUS_STEP_SIZE;
  op[lastop].op.torus.constants = NULL;
  op[lastop].op.torus.cs = NULL;

  ++lastop;
}

| TORUS NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER {
  ++lineno;
  op[lastop].opcode = TORUS;
  op[lastop].op.torus.d[0] = $2;
  op[lastop].op.torus.d[1] = $3;
  op[lastop].op.torus.d[2] = $4;
  op[lastop].op.torus.d[3] = 0;
  op[lastop].op.torus.circle_radius = $5;
  op[lastop].op.torus.torus_radius = $6;
  op[lastop].op.torus.step_size = $7;
  op[lastop].op.torus.constants = NULL;
  op[lastop].op.torus.cs = NULL;

  ++lastop;
}

| TORUS NUMBER NUMBER NUMBER NUMBER NUMBER STRING {
  ++lineno;
  op[lastop].opcode = TORUS;
  op[lastop].op.torus.d[0] = $2;
  op[lastop].op.torus.d[1] = $3;
  op[lastop].op.torus.d[2] = $4;
  op[lastop].op.torus.d[3] = 0;
  op[lastop].op.torus.circle_radius = $5;
  op[lastop].op.torus.torus_radius = $6;
  op[lastop].op.torus.step_size = DEFAULT_TORUS_STEP_SIZE;
  op[lastop].op.torus.constants = NULL;
  m = new_matrix(4,4);
  op[lastop].op.torus.cs = add_symbol($7,SYM_MATRIX,m);
  ++lastop;
}

| TORUS NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER STRING {
  ++lineno;
  op[lastop].opcode = TORUS;
  op[lastop].op.torus.d[0] = $2;
  op[lastop].op.torus.d[1] = $3;
  op[lastop].op.torus.d[2] = $4;
  op[lastop].op.torus.d[3] = 0;
  op[lastop].op.torus.circle_radius = $5;
  op[lastop].op.torus.torus_radius = $6;
  op[lastop].op.torus.step_size = $7;
  op[lastop].op.torus.constants = NULL;
  m = new_matrix(4,4);
  op[lastop].op.torus.cs = add_symbol($8,SYM_MATRIX,m);
  ++lastop;
}

| TORUS STRING NUMBER NUMBER NUMBER NUMBER NUMBER {
  ++lineno;
  op[lastop].opcode = TORUS;
  op[lastop].op.torus.d[0] = $3;
  op[lastop].op.torus.d[1] = $4;
  op[lastop].op.torus.d[2] = $5;
  op[lastop].op.torus.d[3] = 0;
  op[lastop].op.torus.circle_radius = $6;
  op[lastop].op.torus.torus_radius = $7;
  op[lastop].op.torus.step_size = DEFAULT_TORUS_STEP_SIZE;
  op[lastop].op.torus.cs = NULL;
  c = (struct constants *)malloc(sizeof(struct constants));
  op[lastop].op.torus.constants = add_symbol($2,SYM_CONSTANTS,c);

  ++lastop;
}

| TORUS STRING NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER {
  ++lineno;
  op[lastop].opcode = TORUS;
  op[lastop].op.torus.d[0] = $3;
  op[lastop].op.torus.d[1] = $4;
  op[lastop].op.torus.d[2] = $5;
  op[lastop].op.torus.d[3] = 0;
  op[lastop].op.torus.circle_radius = $6;
  op[lastop].op.torus.torus_radius = $7;
  op[lastop].op.torus.step_size = $8;
  op[lastop].op.torus.cs = NULL;
  c = (struct constants *)malloc(sizeof(struct constants));
  op[lastop].op.torus.constants = add_symbol($2,SYM_CONSTANTS,c);

  ++lastop;
}

| TORUS STRING NUMBER NUMBER NUMBER NUMBER NUMBER STRING {
  ++lineno;
  op[lastop].opcode = TORUS;
  op[lastop].op.torus.d[0] = $3;
  op[lastop].op.torus.d[1] = $4;
  op[lastop].op.torus.d[2] = $5;
  op[lastop].op.torus.d[3] = 0;
  op[lastop].op.torus.circle_radius = $6;
  op[lastop].op.torus.torus_radius = $7;
  op[lastop].op.torus.step_size = DEFAULT_TORUS_STEP_SIZE;
  c = (struct constants *)malloc(sizeof(struct constants));
  op[lastop].op.torus.constants = add_symbol($2,SYM_CONSTANTS,c);
  m = new_matrix(4,4);
  op[lastop].op.torus.cs = add_symbol($8,SYM_MATRIX,m);

  ++lastop;
}

| TORUS STRING NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER STRING {
  ++lineno;
  op[lastop].opcode = TORUS;
  op[lastop].op.torus.d[0] = $3;
  op[lastop].op.torus.d[1] = $4;
  op[lastop].op.torus.d[2] = $5;
  op[lastop].op.torus.d[3] = 0;
  op[lastop].op.torus.circle_radius = $6;
  op[lastop].op.torus.torus_radius = $7;
  op[lastop].op.torus.step_size = $8;
  c = (struct constants *)malloc(sizeof(struct constants));
  op[lastop].op.torus.constants = add_symbol($2,SYM_CONSTANTS,c);
  m = new_matrix(4,4);
  op[lastop].op.torus.cs = add_symbol($9,SYM_MATRIX,m);

  ++lastop;
}

| BOX NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER {
  ++lineno;
  op[lastop].opcode = BOX;
  op[lastop].op.box.d0[0] = $2;
  op[lastop].op.box.d0[1] = $3;
  op[lastop].op.box.d0[2] = $4;
  op[lastop].op.box.d0[3] = 0;
  op[lastop].op.box.d1[0] = $5;
  op[lastop].op.box.d1[1] = $6;
  op[lastop].op.box.d1[2] = $7;
  op[lastop].op.box.d1[3] = 0;

  op[lastop].op.box.constants = NULL;
  op[lastop].op.box.cs = NULL;
  ++lastop;
}

| BOX NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER STRING {
  ++lineno;
  op[lastop].opcode = BOX;
  op[lastop].op.box.d0[0] = $2;
  op[lastop].op.box.d0[1] = $3;
  op[lastop].op.box.d0[2] = $4;
  op[lastop].op.box.d0[3] = 0;
  op[lastop].op.box.d1[0] = $5;
  op[lastop].op.box.d1[1] = $6;
  op[lastop].op.box.d1[2] = $7;
  op[lastop].op.box.d1[3] = 0;

  op[lastop].op.box.constants = NULL;
  m = new_matrix(4,4);
  op[lastop].op.box.cs = add_symbol($8,SYM_MATRIX,m);
  ++lastop;
}

| BOX STRING NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER {
  ++lineno;
  op[lastop].opcode = BOX;
  op[lastop].op.box.d0[0] = $3;
  op[lastop].op.box.d0[1] = $4;
  op[lastop].op.box.d0[2] = $5;
  op[lastop].op.box.d0[3] = 0;
  op[lastop].op.box.d1[0] = $6;
  op[lastop].op.box.d1[1] = $7;
  op[lastop].op.box.d1[2] = $8;
  op[lastop].op.box.d1[3] = 0;
  c = (struct constants *)malloc(sizeof(struct constants));
  op[lastop].op.box.constants = add_symbol($2,SYM_CONSTANTS,c);
  op[lastop].op.box.cs = NULL;
  ++lastop;
}

| BOX STRING NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER STRING {
  ++lineno;
  op[lastop].opcode = BOX;
  op[lastop].op.box.d0[0] = $3;
  op[lastop].op.box.d0[1] = $4;
  op[lastop].op.box.d0[2] = $5;
  op[lastop].op.box.d0[3] = 0;
  op[lastop].op.box.d1[0] = $6;
  op[lastop].op.box.d1[1] = $7;
  op[lastop].op.box.d1[2] = $8;
  op[lastop].op.box.d1[3] = 0;
  c = (struct constants *)malloc(sizeof(struct constants));
  op[lastop].op.box.constants = add_symbol($2,SYM_CONSTANTS,c);
  m = new_matrix(4,4);
  op[lastop].op.box.cs = add_symbol($9,SYM_MATRIX,m);

  ++lastop;
}

| LINE NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER {
  ++lineno;
  op[lastop].opcode = LINE;
  op[lastop].op.line.p0[0] = $2;
  op[lastop].op.line.p0[1] = $3;
  op[lastop].op.line.p0[2] = $4;
  op[lastop].op.line.p0[3] = 0;
  op[lastop].op.line.p1[0] = $5;
  op[lastop].op.line.p1[1] = $6;
  op[lastop].op.line.p1[2] = $7;
  op[lastop].op.line.p1[3] = 0;
  op[lastop].op.line.constants = NULL;
  op[lastop].op.line.cs0 = NULL;
  op[lastop].op.line.cs1 = NULL;
  ++lastop;
}

| LINE NUMBER NUMBER NUMBER STRING NUMBER NUMBER NUMBER {
  ++lineno;
  op[lastop].opcode = LINE;
  op[lastop].op.line.p0[0] = $2;
  op[lastop].op.line.p0[1] = $3;
  op[lastop].op.line.p0[2] = $4;
  op[lastop].op.line.p0[3] = 0;
  op[lastop].op.line.p1[0] = $6;
  op[lastop].op.line.p1[1] = $7;
  op[lastop].op.line.p1[2] = $8;
  op[lastop].op.line.p1[3] = 0;
  op[lastop].op.line.constants = NULL;
  m = new_matrix(4,4);
  op[lastop].op.line.cs0 = add_symbol($5,SYM_MATRIX,m);
  op[lastop].op.line.cs1 = NULL;
  ++lastop;
}

| LINE NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER STRING {
  ++lineno;
  op[lastop].opcode = LINE;
  op[lastop].op.line.p0[0] = $2;
  op[lastop].op.line.p0[1] = $3;
  op[lastop].op.line.p0[2] = $4;
  op[lastop].op.line.p0[3] = 0;
  op[lastop].op.line.p1[0] = $5;
  op[lastop].op.line.p1[1] = $6;
  op[lastop].op.line.p1[2] = $7;
  op[lastop].op.line.p1[3] = 0;
  op[lastop].op.line.constants = NULL;
  op[lastop].op.line.cs0 = NULL;
  m = new_matrix(4,4);
  op[lastop].op.line.cs1 = add_symbol($8,SYM_MATRIX,m);
  ++lastop;
}

| LINE NUMBER NUMBER NUMBER STRING NUMBER NUMBER NUMBER STRING {
  ++lineno;
  op[lastop].opcode = LINE;
  op[lastop].op.line.p0[0] = $2;
  op[lastop].op.line.p0[1] = $3;
  op[lastop].op.line.p0[2] = $4;
  op[lastop].op.line.p0[3] = 0;
  op[lastop].op.line.p1[0] = $6;
  op[lastop].op.line.p1[1] = $7;
  op[lastop].op.line.p1[2] = $8;
  op[lastop].op.line.p1[3] = 0;
  op[lastop].op.line.constants = NULL;
  m = new_matrix(4,4);
  op[lastop].op.line.cs0 = add_symbol($5,SYM_MATRIX,m);
  m = new_matrix(4,4);
  op[lastop].op.line.cs1 = add_symbol($9,SYM_MATRIX,m);
  ++lastop;
}

/* now do constants, and constants with the cs stuff */
| LINE STRING NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER {
  ++lineno;
  op[lastop].opcode = LINE;
  op[lastop].op.line.p0[0] = $3;
  op[lastop].op.line.p0[1] = $4;
  op[lastop].op.line.p0[2] = $5;
  op[lastop].op.line.p0[3] = 0;
  op[lastop].op.line.p1[0] = $6;
  op[lastop].op.line.p1[1] = $7;
  op[lastop].op.line.p1[2] = $8;
  op[lastop].op.line.p1[3] = 0;
  c = (struct constants *)malloc(sizeof(struct constants));
  op[lastop].op.line.constants = add_symbol($2,SYM_CONSTANTS,c);
  op[lastop].op.line.cs0 = NULL;
  op[lastop].op.line.cs1 = NULL;
  ++lastop;
}

| LINE STRING NUMBER NUMBER NUMBER STRING NUMBER NUMBER NUMBER {
  ++lineno;
  op[lastop].opcode = LINE;
  op[lastop].op.line.p0[0] = $3;
  op[lastop].op.line.p0[1] = $4;
  op[lastop].op.line.p0[2] = $5;
  op[lastop].op.line.p0[3] = 0;
  op[lastop].op.line.p1[0] = $7;
  op[lastop].op.line.p1[1] = $8;
  op[lastop].op.line.p1[2] = $9;
  op[lastop].op.line.p1[3] = 0;
  c = (struct constants *)malloc(sizeof(struct constants));
  op[lastop].op.line.constants = add_symbol($2,SYM_CONSTANTS,c);
  m = new_matrix(4,4);
  op[lastop].op.line.cs0 = add_symbol($6,SYM_MATRIX,m);
  op[lastop].op.line.cs1 = NULL;
  ++lastop;
}

| LINE STRING NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER STRING {
  ++lineno;
  op[lastop].opcode = LINE;
  op[lastop].op.line.p0[0] = $3;
  op[lastop].op.line.p0[1] = $4;
  op[lastop].op.line.p0[2] = $5;
  op[lastop].op.line.p0[3] = 0;
  op[lastop].op.line.p1[0] = $6;
  op[lastop].op.line.p1[1] = $7;
  op[lastop].op.line.p1[2] = $8;
  op[lastop].op.line.p1[3] = 0;
  c = (struct constants *)malloc(sizeof(struct constants));
  op[lastop].op.line.constants = add_symbol($2,SYM_CONSTANTS,c);
  op[lastop].op.line.cs0 = NULL;
  m = new_matrix(4,4);
  op[lastop].op.line.cs1 = add_symbol($9,SYM_MATRIX,m);
  op[lastop].op.line.cs0 = NULL;
  ++lastop;
}

| LINE STRING NUMBER NUMBER NUMBER STRING NUMBER NUMBER NUMBER STRING {
  ++lineno;
  op[lastop].opcode = LINE;
  op[lastop].op.line.p0[0] = $3;
  op[lastop].op.line.p0[1] = $4;
  op[lastop].op.line.p0[2] = $5;
  op[lastop].op.line.p0[3] = 0;
  op[lastop].op.line.p1[0] = $7;
  op[lastop].op.line.p1[1] = $8;
  op[lastop].op.line.p1[2] = $9;
  op[lastop].op.line.p1[3] = 0;
  c = (struct constants *)malloc(sizeof(struct constants));
  op[lastop].op.line.constants = add_symbol($2,SYM_CONSTANTS,c);
  m = new_matrix(4,4);
  op[lastop].op.line.cs0 = add_symbol($6,SYM_MATRIX,m);
  m = new_matrix(4,4);
  op[lastop].op.line.cs1 = add_symbol($10,SYM_MATRIX,m);
  ++lastop;
}

| MESH CO STRING {
  ++lineno;
  op[lastop].opcode = MESH;
  strcpy(op[lastop].op.mesh.name,$3);
  op[lastop].op.mesh.constants = NULL;
  op[lastop].op.mesh.cs = NULL;
  ++lastop;
}

| MESH STRING CO STRING { /* name and constants */
  ++lineno;
  op[lastop].opcode = MESH;
  strcpy(op[lastop].op.mesh.name,$4);
  c = (struct constants *)malloc(sizeof(struct constants));
  op[lastop].op.mesh.constants = add_symbol($2,SYM_CONSTANTS,c);
  op[lastop].op.mesh.cs = NULL;
  ++lastop;
}

| MESH STRING CO STRING STRING {
  ++lineno;
  op[lastop].opcode = MESH;
  strcpy(op[lastop].op.mesh.name,$4);
  c = (struct constants *)malloc(sizeof(struct constants));
  op[lastop].op.mesh.constants = add_symbol($2,SYM_CONSTANTS,c);
  m = new_matrix(4,4);
  op[lastop].op.mesh.cs = add_symbol($5,SYM_MATRIX,m);
  ++lastop;
}

| SET STRING NUMBER {
  ++lineno;
  op[lastop].opcode = SET;
  op[lastop].op.set.p = add_symbol($2,SYM_STRING,NULL);
  set_value(op[lastop].op.set.p,$3);
  op[lastop].op.set.val = $3;
  ++lastop;
}

| SCALE NUMBER NUMBER NUMBER STRING {
  ++lineno;
  op[lastop].opcode = SCALE;
  op[lastop].op.scale.d[0] = $2;
  op[lastop].op.scale.d[1] = $3;
  op[lastop].op.scale.d[2] = $4;
  op[lastop].op.scale.d[3] = 0;
  op[lastop].op.scale.p = add_symbol($5,SYM_STRING,NULL);
  ++lastop;
}

| SCALE NUMBER NUMBER NUMBER {
  ++lineno;
  op[lastop].opcode = SCALE;
  op[lastop].op.scale.d[0] = $2;
  op[lastop].op.scale.d[1] = $3;
  op[lastop].op.scale.d[2] = $4;
  op[lastop].op.scale.d[3] = 0;
  op[lastop].op.scale.p = NULL;
  ++lastop;
}

| ROTATE STRING NUMBER STRING {
  ++lineno;
  op[lastop].opcode = ROTATE;
  switch (*$2) {
    case 'x':
    case 'X': 
      op[lastop].op.rotate.axis = X_AXIS;
      break;
    case 'y':
    case 'Y': 
      op[lastop].op.rotate.axis = Y_AXIS;
      break;
    case 'z':
    case 'Z': 
      op[lastop].op.rotate.axis = Z_AXIS;
      break;
    }

  op[lastop].op.rotate.degrees = $3;
  op[lastop].op.rotate.p = add_symbol($4,SYM_STRING,NULL);
  
  ++lastop;
}

| ROTATE STRING NUMBER {
  ++lineno;
  op[lastop].opcode = ROTATE;
  switch (*$2)
    {
    case 'x':
    case 'X': 
      op[lastop].op.rotate.axis = 0;
      break;
    case 'y':
    case 'Y': 
      op[lastop].op.rotate.axis = 1;
      break;
    case 'z':
    case 'Z': 
      op[lastop].op.rotate.axis = 2;
      break;
    }
  op[lastop].op.rotate.degrees = $3;
  op[lastop].op.rotate.p = NULL;
  ++lastop;
}

| BASENAME STRING {
  ++lineno;
  if (!is_animation) {
    print_error("Basename can only be set in an animation script. "
                "The frames command must be the first command in the script!");
    free_table();
    exit(EXIT_FAILURE);
  }
  else if (lineno != 2) {
    print_error("Basename must be the second command in an animation script! ");
    free_table();
    exit(EXIT_FAILURE);
  }
  op[lastop].opcode = BASENAME;
  op[lastop].op.basename.p = add_symbol($2,SYM_STRING,NULL);
  ++lastop;
}

| SAVE_KNOBS STRING {
  ++lineno;
  op[lastop].opcode = SAVE_KNOBS;
  op[lastop].op.save_knobs.p = add_symbol($2,SYM_STRING,NULL);
  ++lastop;
}

| TWEEN NUMBER NUMBER STRING STRING {
  ++lineno;
  op[lastop].opcode = TWEEN;
  op[lastop].op.tween.start_frame = $2;
  op[lastop].op.tween.end_frame = $3;
  op[lastop].op.tween.knob_list0 = add_symbol($4,SYM_STRING,NULL);
  op[lastop].op.tween.knob_list1 = add_symbol($5,SYM_STRING,NULL);
  ++lastop;
}

| FRAMES NUMBER {
  ++lineno;
  if (lineno == 1) {
    is_animation = TRUE;
  }
  else {
    print_error("Frames must be the first command in an animation script!");
    free_table();
    exit(EXIT_FAILURE);
  }
  op[lastop].opcode = FRAMES;
  op[lastop].op.frames.num_frames = (int) $2;
  ++lastop;
}

| VARY STRING NUMBER NUMBER NUMBER NUMBER {
  if (!is_animation) {
    print_error("Vary can only be used in an animation script. "
                "The frames command must be the first command in the script!");
    free_table();
    exit(EXIT_FAILURE);
  }
  ++lineno;
  op[lastop].opcode = VARY;
  op[lastop].op.vary.p = add_symbol($2,SYM_STRING,NULL);
  op[lastop].op.vary.start_frame = (int) $3;
  op[lastop].op.vary.end_frame = (int) $4;
  op[lastop].op.vary.start_val = $5;
  op[lastop].op.vary.end_val = $6;
  ++lastop;
}

| PUSH {
  ++lineno;
  op[lastop].opcode = PUSH;
  ++lastop;
}

| GENERATE_RAYFILES {
  ++lineno;
  op[lastop].opcode = GENERATE_RAYFILES;
  ++lastop;
}

| POP {
  ++lineno;
  op[lastop].opcode = POP;
  ++lastop;
}

| SAVE STRING {
  ++lineno;
  op[lastop].opcode = SAVE;
  op[lastop].op.save.p = add_symbol($2,SYM_STRING,NULL);
  ++lastop;
}

| SETKNOBS NUMBER {
  ++lineno;
  op[lastop].opcode = SETKNOBS;
  op[lastop].op.setknobs.value = $2;
  ++lastop;
}

| FOCAL NUMBER {
  ++lineno;
  op[lastop].opcode = FOCAL;
  op[lastop].op.focal.value = $2;
  ++lastop;
}

| DISPLAY {
  ++lineno;
  op[lastop].opcode = DISPLAY;
  ++lastop;
}

| WEB {
  ++lineno;
  op[lastop].opcode = WEB;
  ++lastop;
}

| AMBIENT NUMBER NUMBER NUMBER {
  ++lineno;
  op[lastop].opcode = AMBIENT;
  op[lastop].op.ambient.c[0] = $2;
  op[lastop].op.ambient.c[1] = $3;
  op[lastop].op.ambient.c[2] = $4;
  ++lastop;
};

%%

/* ================ END RULES ================ */

/* ================ SUBROUTINES ============== */
void yyerror(const char *s) {
    fprintf(stderr, "Error on line %d: %s\n", lineno, s);
}

int yywrap() {
    return 1;
}

static void sighandler(int signo) {
    if (signo == SIGINT) {
        // TODO free vary_nodes and others
        free_table();
        exit(EXIT_SUCCESS);
    }
}

int main(int argc, char *argv[]) {
    if (argc > 1) {
        yyin = fopen(argv[1], "r");

        if (yyin == NULL) {
            print_error("Unable to open input stream.");
            exit(EXIT_FAILURE);
        }

    }

    int ret = yyparse();
    if (ret == 1) { // Syntax error
        print_error("Fatal error.");
        exit(EXIT_FAILURE);
    }

    #ifdef DEBUG
    print_pcode();
    #endif
    if (is_animation) {
        parse_animation_cmds();
        exec_animation();
    }
    else {
        exec(FALSE);
    }

    free_table();
    return 0;    
}

// TODO dynamically allocated op array size
/* ============= END SUBROUTINES ============= */
