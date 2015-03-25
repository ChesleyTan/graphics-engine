#include "repl.h"

// REPL variables
char keep_alive = 1;
char input[INPUT_BUF_SIZE];
int rl_child_pid;
const char *prompt = ">> ";
int pipes[2];
int history_pipes[2];

screen s;
color c;
struct matrix *mat;
struct matrix *trans_mat;

static void readline_sigint_handler() {
    // Exit gracefully when killed with SIGINT
    exit(SIGINT_EXIT_CODE);
}

static void sighandler(int signo) {
    if (signo == SIGINT) {
        if (rl_child_pid) {
            // Kill readline process to refresh prompt
            kill(rl_child_pid, SIGINT);
        }
    }
}

int main(int argc, char *argv[]) {
    signal(SIGINT, sighandler);
    s = new_screen();
    mat = new_matrix(4, 1);
    trans_mat = new_matrix(4, 4);
    synchronize_variables(NULL, trans_mat, mat, s);
    while(keep_alive) {
        if (pipe(pipes) < 0) {
            print_error("Could not open pipe for REPL.");
            exit(1);
        }
        if (pipe(history_pipes) < 0) {
            print_error("Could not open pipe for command history.");
            exit(1);
        }
        // Fork process for readline
        rl_child_pid = fork();
        if (!rl_child_pid) {
            signal(SIGINT, readline_sigint_handler);
            close(pipes[0]);
            close(history_pipes[0]);
            char *line = readline(prompt);
            if (line == NULL) {
                printf("\n[Reached EOF]\n");
                // Free dynamically allocated memory before exiting
                free(line);
                // Close pipes before exiting
                close(pipes[1]);
                close(history_pipes[1]);
                exit(EOF_EXIT_CODE);
            }
            size_t read_size = strlen(line) + 1; // Include null terminator
            // Limit write size to INPUT_BUF_SIZE
            size_t write_size = (INPUT_BUF_SIZE > read_size)
                                ? read_size : INPUT_BUF_SIZE;
            line[write_size - 1] = '\0';
            // Write input to pipe -> parser
            write(pipes[1], line, write_size);
            // Write input to pipe -> readline history
            write(history_pipes[1], line, write_size);
            // Free dynamically allocated memory before exiting
            free(line);
            // Close pipes before exiting
            close(pipes[1]);
            close(history_pipes[1]);
            exit(0);
        }
        else {
            int status;
            waitpid(rl_child_pid, &status, 0);
            if (WIFEXITED(status)) {
                status = WEXITSTATUS(status);
                // Exit if EOF reached
                if (status == EOF_EXIT_CODE) {
                    // Close pipes before exiting
                    close(pipes[0]);
                    close(pipes[1]);
                    close(history_pipes[0]);
                    close(history_pipes[1]);
                    // Clear readline's history
                    clear_history();
                    // Free dynamically allocated memory
                    free_variables();
                    exit(0);
                }
                // Prepare for re-displaying prompt after Ctrl-c
                else if (status == SIGINT_EXIT_CODE) {
                    // Close pipes before relooping
                    close(pipes[0]);
                    close(pipes[1]);
                    close(history_pipes[0]);
                    close(history_pipes[1]);
                    // Go to new line before relooping
                    write(STDOUT_FILENO, "\n", 1);
                    continue;
                }
            }
            close(pipes[1]);
            close(history_pipes[1]);
            // Read user input for adding to history
            int bytes = read(history_pipes[0], input, INPUT_BUF_SIZE);
            close(history_pipes[0]);
            // Add the input to readline's history
            add_history(input);
            synchronize_variables(fdopen(pipes[0], "r"), NULL, NULL, NULL);
            parse_file(0);
            close(pipes[0]);
        }
    }
    return 0;
}
