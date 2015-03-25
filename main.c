#include "main.h"

static void sighandler(int signo) {
    if (signo == SIGINT) {
        free_variables();
        exit(0);
    }
}

int main(int argc, char *argv[]) {
    signal(SIGINT, sighandler);
    if (argc > 1) {
        #ifdef DEBUG
        print_debug("Got file argument: %s", argv[1]);
        #endif
        FILE *file = fopen(argv[1], "r");
        if (file == NULL) {
            print_error("Could not open file for reading: %s", argv[1]);
            exit(1);
        }
        struct matrix *mat = new_matrix(4, 1);
        screen s = new_screen();
        struct matrix *trans_mat = new_matrix(4,4);
        synchronize_variables(file, trans_mat, mat, s);
        parse_file(1);
        free_variables();
        exit(0);
    }

    printf("Usage: ./main [script-file]\n");

    return 0;
}
// vim: ts=4:et:sts:sw=4:sr
