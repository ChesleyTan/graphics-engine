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

struct vary_node {
    SYMTAB *knob;
    double next_value;
    struct vary_node *next;
};

// TODO update doc
extern int num_frames;
extern char *basename;
// Array of vary_node linked-lists
extern struct vary_node **vary_knobs;

// TODO update doc
screen exec(char return_screen);
void parse_animation_cmds();
void exec_animation();
struct vary_node *get_vary_knobs_tail(int frame);
int vary_node_uniq(int frame, SYMTAB *knob);
void free_vary_knobs();
void free_all();
#endif
