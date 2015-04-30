#pragma once
#ifndef REPL_H
#define REPL_H
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <signal.h>
#include <sys/wait.h>
#include <readline/readline.h>
#include <readline/history.h>

#include "parser-old.h"

#define INPUT_BUF_SIZE 512
#define EOF_EXIT_CODE 10
#define SIGINT_EXIT_CODE 11

/*======== static void readline_sigint_handler() ===========
Inputs:
Returns:
Handles the graceful exit of the readline child process when SIGINT is received.
==========================================================*/
static void readline_sigint_handler();

/*======== static void readline_sigint_handler() ===========
Inputs:
Returns: 
When SIGINT is received, the readline child process is killed so that
the prompt can be refreshed and the input buffer cleared.
==========================================================*/
static void sighandler(int signo);

#endif
