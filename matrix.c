#include "matrix.h"

struct matrix *new_matrix(int rows, int cols) {
    double **tmp;
    int i;
    struct matrix *m;

    if (rows < 1) {
        rows = 1;
    }
    if (cols < 1) {
        cols = 1;
    }

    tmp = (double **) malloc(rows * sizeof(double *));
    if (tmp == NULL) {
        print_error("Memory allocation error.");
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < rows; ++i) {
        tmp[i] = (double *) calloc(cols, sizeof(double));
        if (tmp[i] == NULL) {
            print_error("Memory allocation error.");
            exit(EXIT_FAILURE);
        }
    }
    m = (struct matrix *) malloc(sizeof(struct matrix));
    if (m == NULL) {
        print_error("Memory allocation error.");
        exit(EXIT_FAILURE);
    }
    m->m = tmp;
    m->rows = rows;
    m->cols = cols;
    m->lastcol = 0;

    return m;
}

void free_matrix(struct matrix *m) {
    if (m != NULL) {
        int i;
        for (i = 0; i < m->rows; ++i) {
            free(m->m[i]);
        }
        free(m->m);
        free(m);
    }
}

void grow_matrix(struct matrix *m, int newcols) {
    int i;
    for (i = 0; i < m->rows; ++i) {
        double *ptr = realloc(m->m[i], newcols * sizeof(double));
        if (ptr != NULL) {
            m->m[i] = ptr;
        }
        else {
            print_error("Memory reallocation error.");
            exit(EXIT_FAILURE);
        }
    }
    m->cols = newcols;
}

void copy_matrix(struct matrix *a, struct matrix *b) {
    int r, c;
    for (r = 0; r < a->rows; ++r) {
        for (c = 0; c < a->cols; ++c) {
            b->m[r][c] = a->m[r][c];
            b->lastcol = a->lastcol;
        }
    }
}

void print_matrix(struct matrix *m) {
    double **pts = m->m;
    int i,u;
    for (u = 0; u < m->lastcol; ++u) {
        printf("(");
        for (i = 0; i < m->rows; ++i) {
            printf("%lf", pts[i][u]);
            if (i < m->rows - 1) {
                printf(",");
            }
        }
        printf(")\n");
    }
}

void ident(struct matrix *m) {
    int r, c;
    for (r = 0; r < m->rows; ++r) {
        for (c = 0; c < m->cols; ++c)  {
            if ( r == c ) {
                m->m[r][c] = 1;
            }
            else {
                m->m[r][c] = 0;
            }
        }
    }
    m->lastcol = m->cols;
}

void scalar_mult(double x, struct matrix *m) {
    int r,c;
    for (r = 0; r < m->rows; ++r) {
        for (c = 0; c < m->lastcol; ++c) {
            m->m[r][c] *= x;
        }
    }
}

struct matrix *matrix_mult(struct matrix *a, struct matrix *b) {
    if (a->lastcol != b->rows) {
        print_error("Column-row number mismatch in matrix multiplication: "
                    "(%d x %d) * (%d x %d)",
                    a->rows, a->lastcol, b->rows, b->lastcol);
        exit(EXIT_FAILURE);
    }
    struct matrix *m = new_matrix(a->rows, b->cols);
    m->lastcol = b->lastcol;
    int i,j,k;
    double sum;
    for (i = 0; i < a->rows; ++i) {
        for (j = 0; j < b->lastcol; ++j) {
            sum = 0;
            for (k = 0; k < b->rows; ++k) {
                sum += a->m[i][k] * b->m[k][j];
            }
            m->m[i][j] = sum;
        }
    }
    return m;
}

struct matrix *make_translate(double x, double y, double z) {
    /*
    |1 0 0 a||x|   |x + a|
    |0 1 0 b||y|   |y + b|
    |0 0 1 c||z| = |z + c|
    |0 0 0 1||1|   |  1  |
    */
    struct matrix *m = new_matrix(4, 4);
    m->m[0][0] = 1;
    m->m[0][1] = 0;
    m->m[0][2] = 0;
    m->m[0][3] = x;
    m->m[1][0] = 0;
    m->m[1][1] = 1;
    m->m[1][2] = 0;
    m->m[1][3] = y;
    m->m[2][0] = 0;
    m->m[2][1] = 0;
    m->m[2][2] = 1;
    m->m[2][3] = z;
    m->m[3][0] = 0;
    m->m[3][1] = 0;
    m->m[3][2] = 0;
    m->m[3][3] = 1;
    m->lastcol = 4;
    return m;
}

struct matrix *make_scale(double x, double y, double z) {
    /*
    |a 0 0 0||x|   |ax|
    |0 b 0 0||y|   |by|
    |0 0 c 0||z| = |cz|
    |0 0 0 1||1|   |1 |
    */
    struct matrix *m = new_matrix(4, 4);
    m->m[0][0] = x;
    m->m[0][1] = 0;
    m->m[0][2] = 0;
    m->m[0][3] = 0;
    m->m[1][0] = 0;
    m->m[1][1] = y;
    m->m[1][2] = 0;
    m->m[1][3] = 0;
    m->m[2][0] = 0;
    m->m[2][1] = 0;
    m->m[2][2] = z;
    m->m[2][3] = 0;
    m->m[3][0] = 0;
    m->m[3][1] = 0;
    m->m[3][2] = 0;
    m->m[3][3] = 1;
    m->lastcol = 4;
    return m;
}

struct matrix *make_rotX(double theta) {
    /*
    |1       0        0       0||x|   |            x            |
    |0 cos(theta) -sin(theta) 0||y|   |ycos(theta) - zsin(theta)|
    |0 sin(theta)  cos(theta) 0||z| = |ysin(theta) + zcos(theta)|
    |0       0        0       1||1|   |            1            |
    */
    struct matrix *m = new_matrix(4, 4);
    theta = deg_to_rad(theta);
    m->m[0][0] = 1;
    m->m[0][1] = 0;
    m->m[0][2] = 0;
    m->m[0][3] = 0;
    m->m[1][0] = 0;
    m->m[1][1] = cos(theta);
    m->m[1][2] = -1*sin(theta);
    m->m[1][3] = 0;
    m->m[2][0] = 0;
    m->m[2][1] = sin(theta);
    m->m[2][2] = cos(theta);
    m->m[2][3] = 0;
    m->m[3][0] = 0;
    m->m[3][1] = 0;
    m->m[3][2] = 0;
    m->m[3][3] = 1;
    m->lastcol = 4;
    return m;
}

struct matrix *make_rotY(double theta) {
    /*
    |cos(theta) 0 -sin(theta) 0||x|   |xcos(theta) - zsin(theta)|
    |0       1        0       0||y|   |            y            |
    |sin(theta) 0  cos(theta) 0||z| = |xsin(theta) + zcos(theta)|
    |0       0        0       1||1|   |            1            |
    */
    struct matrix *m = new_matrix(4, 4);
    theta = deg_to_rad(theta);
    m->m[0][0] = cos(theta);
    m->m[0][1] = 0;
    m->m[0][2] = -1*sin(theta);
    m->m[0][3] = 0;
    m->m[1][0] = 0;
    m->m[1][1] = 1;
    m->m[1][2] = 0;
    m->m[1][3] = 0;
    m->m[2][0] = sin(theta);
    m->m[2][1] = 0;
    m->m[2][2] = cos(theta);
    m->m[2][3] = 0;
    m->m[3][0] = 0;
    m->m[3][1] = 0;
    m->m[3][2] = 0;
    m->m[3][3] = 1;
    m->lastcol = 4;
    return m;
}

struct matrix *make_rotZ(double theta) {
    /*
    |cos(theta) -sin(theta) 0 0||x|   |xcos(theta) - ysin(theta)|
    |sin(theta) cos(theta)  0 0||y|   |ycos(theta) + xsin(theta)|
    |0       0        1       0||z| = |            z            |
    |0       0        0       1||1|   |            1            |
    */
    struct matrix *m = new_matrix(4, 4);
    theta = deg_to_rad(theta);
    m->m[0][0] = cos(theta);
    m->m[0][1] = -1*sin(theta);
    m->m[0][2] = 0;
    m->m[0][3] = 0;
    m->m[1][0] = sin(theta);
    m->m[1][1] = cos(theta);
    m->m[1][2] = 0;
    m->m[1][3] = 0;
    m->m[2][0] = 0;
    m->m[2][1] = 0;
    m->m[2][2] = 1;
    m->m[2][3] = 0;
    m->m[3][0] = 0;
    m->m[3][1] = 0;
    m->m[3][2] = 0;
    m->m[3][3] = 1;
    m->lastcol = 4;
    return m;
}

struct matrix *apply_transform(struct matrix *trans_mat,
                               struct matrix *points) {
    return matrix_mult(trans_mat, points);
}

struct matrix *make_hermite_coefficients() {
    /*
    |  2 -2  1  1 || P0 |   | a |
    | -3  3 -2 -1 || P1 | = | b |
    |  0  0  1  0 || R0 |   | c |
    |  1  0  0  0 || R1 |   | d |
    */
    struct matrix *coeff = new_matrix(4, 4);
    coeff->m[0][0] = 2;
    coeff->m[0][1] = -2;
    coeff->m[0][2] = 1;
    coeff->m[0][3] = 1;
    coeff->m[1][0] = -3;
    coeff->m[1][1] = 3;
    coeff->m[1][2] = -2;
    coeff->m[1][3] = -1;
    coeff->m[2][0] = 0;
    coeff->m[2][1] = 0;
    coeff->m[2][2] = 1;
    coeff->m[2][3] = 0;
    coeff->m[3][0] = 1;
    coeff->m[3][1] = 0;
    coeff->m[3][2] = 0;
    coeff->m[3][3] = 0;
    coeff->lastcol = 4;
    return coeff;
}

struct matrix *make_bezier_coefficients() {
    /*
    | -1  3 -3 1 |   | P0 |   | a |
    |  3 -6  3 0 | * | P1 | = | b |
    | -3  3  0 0 |   | P2 |   | c |
    |  1  0  0 0 |   | P3 |   | d |
    */
    struct matrix *coeff = new_matrix(4, 4);
    coeff->m[0][0] = -1;
    coeff->m[0][1] = 3;
    coeff->m[0][2] = -3;
    coeff->m[0][3] = 1;
    coeff->m[1][0] = 3;
    coeff->m[1][1] = -6;
    coeff->m[1][2] = 3;
    coeff->m[1][3] = 0;
    coeff->m[2][0] = -3;
    coeff->m[2][1] = 3;
    coeff->m[2][2] = 0;
    coeff->m[2][3] = 0;
    coeff->m[3][0] = 1;
    coeff->m[3][1] = 0;
    coeff->m[3][2] = 0;
    coeff->m[3][3] = 0;
    coeff->lastcol = 4;
    return coeff;

}
// vim: ts=4:et:sts:sw=4:sr
