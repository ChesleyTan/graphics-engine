#pragma once
#ifndef STACK_H
#define STACK_H

#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

#define STACK_SIZE 2

struct stack {
    struct matrix **data;
    int size;
    int top;
};

struct stack * new_stack();
void push(struct stack *s);
void pop(struct stack *s);

void free_stack(struct stack *s);
void print_stack(struct stack *s);

#endif
