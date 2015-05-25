#pragma once
#ifndef UTILS_H
#define UTILS_H

#define FALSE 0
#define TRUE 1

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <errno.h>
#include <string.h>
#include <math.h>


// ANSI Escape codes for output formatting
static const char *reset = "\e[0m";
static const char *bold_prefix = "\e[1;";
static const char *dim_prefix = "\e[2;";
static const char *underline_prefix = "\e[4;";
static const char *prefix = "\e[";
static const char *fg_red_160 = "38;5;160m";
static const char *fg_red_196 = "38;5;196m";
static const char *fg_yellow_220 = "38;5;220m";
static const char *fg_green_34 = "38;5;34m";
static const char *fg_green_118 = "38;5;118m";
static const char *fg_blue_24 = "38;5;24m";
static const char *fg_blue_30 = "38;5;30m";
static const char *fg_blue_39 = "38;5;39m";
static const char *fg_white = "38;5;15m";
static const char *bg_black = "48;5;0m";

/*======== void print_error() ===========
Inputs:     const char *format,
            ...
Returns:
Prints a red error message to stderr using a variable number of arguments and
the format parameter.
========================================*/
void print_error(const char *format, ...);

/*======== void print_warning() ===========
Inputs:     const char *format,
            ...
Returns:
Prints a yellow warning message to stderr using a variable number of arguments and
the format parameter.
========================================*/
void print_warning(const char *format, ...);

/*======== void print_debug() ===========
Inputs:     const char *format,
            ...
Returns:
Prints a yellow debug message to stderr using a variable number of arguments and
the format parameter.
========================================*/
void print_debug(const char *format, ...);

/*======== void print_errno() ===========
Inputs:     const char *message
Returns:
Prints a red error message to stderr using errno and its associated error
message in addition to the provided message parameter.
========================================*/
void print_errno(const char *message);

/*======== double deg_to_rad() ===========
Inputs:     double degrees
Returns:
The radian value of an angle of "degrees" degrees.
================================================*/
double deg_to_rad(double degrees);

/*======== int *remove_trailing_newline() ===========
Inputs:     char *str
Returns:
Truncates the string at its first newline character.
Returns 0 if a newline was found and replaced with a null character.
Returns -1 if no newline was found.
================================================*/
int remove_trailing_newline(char *str);
#endif
// vim: ts=4:et:sts:sw=4:sr
