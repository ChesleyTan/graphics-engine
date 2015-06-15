#pragma once
#ifndef EXEC_H
#define EXEC_H

#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <unistd.h>
#include "parser.h"
#include "display.h"
#include "draw.h"
#include "stack.h"
#include "symtab.h"
#include "y.tab.h"

/*
 * vary_node holds a pointer to a knob in the SYMTAB, the next value for the
 * knob, and a pointer to the next vary_node in its linked list
 */
struct vary_node {
    SYMTAB *knob;
    double next_value;
    struct vary_node *next;
};

/* _screen holds all the colors/pixels that have been drawn */
extern screen _screen;
/* num_frames holds the number of frames specified for the current animation */
extern int num_frames;
/* basename specifies the base name of the frame image files */
extern char *basename;
/* Array of vary_node linked-lists */
extern struct vary_node **vary_knobs;

/*======== screen exec() ==========
Inputs:     char return_screen
Returns:
Executes the mdl commands (opcodes in the global op array).
If return_screen is TRUE (i.e. not 0), then the screen that results from
executing the mdl commands is returned.
=================================*/
screen exec(char return_screen);

/*======== void parse_animation_cmds() ==========
Inputs:
Returns:
Parses the animation-related mdl commands and sets the necessary knobs and other
variables needed for later execution.
===============================================*/
void parse_animation_cmds();

/*======== void exec_animation() ==========
Inputs:
Returns:
Creates an animation by iteratively executing the mdl commands and saving the
resulting frames to a `frames/` folder in the current directory.
=========================================*/
void exec_animation();

/*======== void get_vary_knobs_tail() ==========
Inputs:     int frame
Returns:
A utility function to get the tail of the `vary_knobs` array linked list for the
specified frame.
==============================================*/
struct vary_node *get_vary_knobs_tail(int frame);

/*======== char vary_node_uniq() ==========
Inputs:     int frame
            SYMTAB *knob
Returns:
A utility function to determine if a SYMTAB *knob is not already in the
`vary_node` linked list for the given frame number. Returns FALSE (0) if the
knob is not in a vary_node, and TRUE (1) otherwise.
=========================================*/
char vary_node_uniq(int frame, SYMTAB *knob);

/*======== char allocate_vary_knobs() ==========
Inputs:
Returns:
A utility function to allocate and initialize the `vary_knobs` array.
==============================================*/
void allocate_vary_knobs();

/*======== char free_exec_screen() ==========
Inputs:
Returns:
A utility function to free the `_screen` used in exec().
===========================================*/
void free_exec_screen();

/*======== char free_vary_knobs() ==========
Inputs:
Returns:
A utility function to free the `vary_knobs` array.
==========================================*/
void free_vary_knobs();

/*======== char free_all() ==========
Inputs:
Returns:
A utility function to free all known currently allocated memory. This should be
called prior to exiting.
===================================*/
void free_all();
#endif
