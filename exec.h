#pragma once
#ifndef EXEC_H
#define EXEC_H
#include "parser.h"
#include "display.h"
#include "draw.h"
#include "stack.h"
#include "y.tab.h"

void exec();
static void draw(screen s, struct matrix *pts, color c);
#endif
