%option stack
%{
#include "y.tab.h"
#include "utils.h"
%}
    /* =============== DEFINITIONS ============= */

    /* Parse state for comments */
%x COMMENT_STATE

ID         [a-zA-Z][a-zA-Z0-9_]*
WHITESPACE [ \t\n]
    /* ============= END DEFINITIONS ============= */

    /* ================== RULES ================== */
%%

\-?[0-9]+ |
\-?[0-9]+\. |
\-?[0-9]+\.[0-9]+ |
\-?\.[0-9]+ {
    yylval.val = atof(yytext);
    return NUMBER;
}

    /* ============= Ignore comments ============= */
"//".* { return COMMENT; }
"/*"                                { yy_push_state(COMMENT_STATE); }
    /* Non-greedy regex */
<COMMENT_STATE>(("*"[^/])?|[^*])*   ;
<COMMENT_STATE>"*/"                 {
                                      yy_pop_state();
                                      return COMMENT;
                                    }
    /* =========================================== */

"draw-mode" { return DRAW_MODE; }
"lines"|"polygons" {
    strncpy(yylval.string, yytext, sizeof(yylval.string) - 1);
    yylval.string[sizeof(yylval.string) - 1] = '\0';
    return DRAW_TYPE;
}
"render-mode" { return RENDER_MODE; }
"wireframe"|"surface" {
    strncpy(yylval.string, yytext, sizeof(yylval.string) - 1);
    yylval.string[sizeof(yylval.string) - 1] = '\0';
    return RENDER_TYPE;
}
"resize" { return RESIZE; }
"light" { return LIGHT; }
"constants" { return CONSTANTS; }
"save_coord_system" { return SAVE_COORDS; }
"camera" { return CAMERA; }
"ambient" { return AMBIENT; }

"torus" { return TORUS; }
"sphere" { return SPHERE; }
"box" { return BOX; }
"line" { return LINE; }
"mesh" { return MESH; }
"texture" { return TEXTURE; }

"set" { return SET; }
"move" { return MOVE; }
"scale" { return SCALE; }
"rotate" { return ROTATE; }
"basename" { return BASENAME; }
"save_knobs" { return SAVE_KNOBS; }
"tween" { return TWEEN; }
"frames" { return FRAMES; }
"vary" { return VARY; }

"push" { return PUSH; }
"pop" { return POP; }
"save" { return SAVE; }
"generate_rayfiles" { return GENERATE_RAYFILES; }

"shade-mode" { return SHADE_MODE; }

"flat"|"gouraud"|"phong"|"raytrace" {
    strncpy(yylval.string, yytext, sizeof(yylval.string) - 1);
    yylval.string[sizeof(yylval.string) - 1] = '\0';
    return SHADING_TYPE;
}

"setknobs" { return SETKNOBS; }
"focal" { return FOCAL; }
"display" { return DISPLAY; }
"web" { return WEB; }

":" { return CO; }

[a-zA-Z][-\.a-zA-Z0-9_]* {
    strncpy(yylval.string, yytext, sizeof(yylval.string) - 1);
    yylval.string[sizeof(yylval.string) - 1] = '\0';
    return STRING;
}

    /* ============= Ignore whitespace =========== */
{WHITESPACE}           ;
    /* =========================================== */
    /* ====== Throw error for anything else ====== */
.   { print_error("Invalid character %s", yytext); }
    /* =========================================== */


%%
    /* ================ END RULES ================ */

