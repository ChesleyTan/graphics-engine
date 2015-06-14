#include "display.h"

double **z_buffer;

screen resize_screen(screen s, int x_res, int y_res) {
    free_screen(s);
    free_z_buffer();
    XRES = x_res;
    YRES = y_res;
    XRES_CARTESIAN = x_res / 2;
    YRES_CARTESIAN = y_res / 2;
    allocate_z_buffer();
    return new_screen();
}

screen new_screen() {
    void *ptr;
    screen s;
    int i;
    ptr = malloc(sizeof(struct point_t *) * XRES);
    if (ptr != NULL) {
        s = (struct point_t **) ptr;
    }
    else {
        print_error("Memory allocation error.");
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < XRES; ++i) {
        ptr = calloc(YRES, sizeof(struct point_t));
        if (ptr != NULL) {
            s[i] = (struct point_t *) ptr;
        }
        else {
            print_error("Memory allocation error.");
            exit(EXIT_FAILURE);
        }
    }
    return s;
}

void clear_screen(screen s) {
    int x, y;
    color c;

    c.red = 0;
    c.green = 0;
    c.blue = 0;

    for (x = 0; x < XRES; ++x) {
        for (y = 0; y < YRES; ++y) {
            s[x][y] = c;
        }
    }
}

void free_screen(screen s) {
    if (s != NULL) {
        int i;
        for (i = 0; i < XRES; ++i) {
            free(s[i]);
        }
        free(s);
    }
}

void plot_absolute(screen s, color c, int x, int y, int z) {
    int newY = YRES - 1 - y;
    if (x >= 0 && x < XRES && newY >= 0 && newY < YRES) {
        if (z_buffer[x][newY] <= z) {
            s[x][newY] = c;
            z_buffer[x][newY] = z;
        }
    }
}

void plot_cartesian(screen s, color c, int x, int y, int z) {
    int newY = YRES_CARTESIAN - y;
    int newX = x + XRES_CARTESIAN;
    if (newX >= 0 && newX < XRES && newY >= 0 && newY < YRES) {
        if (z_buffer[newX][newY] <= z) {
            s[newX][newY] = c;
            z_buffer[newX][newY] = z;
        }
    }
}

void plot_cartesian_wrap(screen s, color c, int x, int y, int z) {
    int newY = YRES_CARTESIAN - y;
    int newX = x + XRES_CARTESIAN;
    if (newY >= YRES) {
        newY %= YRES;
    }
    if (newY < 0) {
        newY = YRES_CARTESIAN - (newY % YRES_CARTESIAN);
    }
    if (newX >= XRES) {
        newX %= XRES;
    }
    if (newX < 0) {
        newX = newX % XRES + XRES_CARTESIAN;
    }
    if (z_buffer[newX][newY] <= z) {
        s[newX][newY] = c;
        z_buffer[newX][newY] = z;
    }
}

void save_ppm(screen s, char *file) {
    int x, y;
    FILE *f;
    color c;

    f = fopen(file, "w");
    if (f == NULL) {
        print_error("Could not open file for writing.");
        exit(EXIT_FAILURE);
    }
    fprintf(f, "P3\n%d %d\n%d\n", XRES, YRES, MAX_COLOR);
    for (y=0; y < YRES; ++y) {
        for (x=0; x < XRES; ++x) {
            c = s[x][y];
            fprintf(f, "%d %d %d ", (int)c.red, (int)c.green, (int)c.blue);
        }
        fprintf(f, "\n");
    }
    fclose(f);
}

void save_extension(screen s, char *file) {
    int x, y;
    FILE *f;
    color c;
    char line[256];

    snprintf(line, sizeof(line), "convert - %s", file);
    line[sizeof(line) - 1] = '\0';

    f = popen(line, "w");
    if (f == NULL) {
        print_error("Could not open file for writing.");
        exit(EXIT_FAILURE);
    }
    fprintf(f, "P3\n%d %d\n%d\n", XRES, YRES, MAX_COLOR);
    for (y=0; y < YRES; ++y) {
        for (x=0; x < XRES; ++x) {
            c = s[x][y];
            fprintf(f, "%d %d %d ", (int)c.red, (int)c.green, (int)c.blue);
        }
        fprintf(f, "\n");
    }
    pclose(f);
}

void display(screen s) {
    int x, y;
    FILE *f;
    color c;

    f = popen("display", "w");
    if (f == NULL) {
        print_error("Could not open file for writing.");
        exit(EXIT_FAILURE);
    }

    fprintf(f, "P3\n%d %d\n%d\n", XRES, YRES, MAX_COLOR);
    for (y=0; y < YRES; ++y) {
        for (x=0; x < XRES; ++x) {
            c = s[x][y];
            fprintf(f, "%d %d %d ", (int)c.red, (int)c.green, (int)c.blue);
        }
        fprintf(f, "\n");
    }
    pclose(f);
}

void allocate_z_buffer() {
    int i, u;
    // Allocate and initialize z_buffer to the minimum double
    z_buffer = (double **) malloc(sizeof(double *) * XRES);
    for (i = 0; i < XRES; ++i) {
        z_buffer[i] = (double *) malloc(sizeof(double) * YRES);
    }
    for (i = 0; i < XRES; ++i) {
        for (u = 0; u < YRES; ++u) {
            z_buffer[i][u] = -DBL_MAX;
        }
    }
}

void free_z_buffer() {
    if (z_buffer != NULL) {
        int i;
        // Free allocated memory for z_buffer
        for (i = 0; i < XRES; ++i) {
            free(z_buffer[i]);
        }
        free(z_buffer);
        z_buffer = NULL;
    }
}

color avg_color(color c1, color c2) {
    color c;
    c.red = (c1.red + c2.red) / 2;
    c.green = (c1.green + c2.green) / 2;
    c.blue = (c1.blue + c2.blue) / 2;
    return c;
}

color add_color(color c1, color c2) {
    color c;
    c.red = c1.red + c2.red;
    c.green = c1.green + c2.green;
    c.blue = c1.blue + c2.blue;
    return c;
}

color subtract_color(color c1, color c2) {
    color c;
    c.red = c1.red - c2.red;
    c.green = c1.green - c2.green;
    c.blue = c1.blue - c2.blue;
    return c;
}

color divide_color(color c, int n) {
    color new_c;
    if (n == 0) {
        new_c.red = 0;
        new_c.green = 0;
        new_c.blue = 0;
    }
    else {
        new_c.red = c.red / n;
        new_c.green = c.green / n;
        new_c.blue = c.blue / n;
    }
    return new_c;
}

color constrain_color(color c) {
    if (c.red > MAX_COLOR) c.red = MAX_COLOR;
    else if (c.red < 0) c.red = 0;
    else if (c.red != c.red) c.red = 0; // NaN check
    if (c.green > MAX_COLOR) c.green = MAX_COLOR;
    else if (c.green < 0) c.green = 0;
    else if (c.green != c.green) c.green = 0;
    if (c.blue > MAX_COLOR) c.blue = MAX_COLOR;
    else if (c.blue < 0) c.blue = 0;
    else if (c.blue != c.blue) c.blue = 0;
    return c;
}

// vim: ts=4:et:sts:sw=4:sr
