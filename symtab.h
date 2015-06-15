#pragma once
#ifndef SYMTAB_H
#define SYMTAB_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "matrix.h"

#define MAX_SYMBOLS 512
#define SYM_MATRIX 1
#define SYM_VALUE 2
#define SYM_CONSTANTS 3
#define SYM_STRING 4

struct constants {
    int iar;
    int iag;
    int iab;
    double kar;
    double kag;
    double kab;
    int idr;
    int idg;
    int idb;
    double kdr;
    double kdg;
    double kdb;
    int isr;
    int isg;
    int isb;
    double ksr;
    double ksg;
    double ksb;
    int spec_expt;
};

typedef struct {
    union{
        struct matrix *m;
        struct constants *c;
        double value;
    } s;
    int type;
    char *name;
} SYMTAB;

extern SYMTAB symtab[MAX_SYMBOLS];
extern int lastsym;

SYMTAB *lookup_symbol(char *name);
SYMTAB *add_symbol(char *name, int type, void *data);
void print_constants(struct constants *p);
void print_symtab();
SYMTAB *add_symbol(char *name, int type, void *data);
void set_value(SYMTAB *p, double value);
void free_table();

#endif
