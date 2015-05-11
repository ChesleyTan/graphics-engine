#include "symtab.h"

// TODO use a hash table
SYMTAB symtab[MAX_SYMBOLS];
int lastsym = 0;

void print_constants(struct constants *p) {
    printf("\tRed -\t  Ka: %6.2f Kd: %6.2f Ks: %6.2f\n",
            p->r[0],p->r[1],p->r[2]);

    printf("\tGreen -\t  Ka: %6.2f Kd: %6.2f Ks: %6.2f\n",
            p->g[0],p->g[1],p->g[2]);

    printf("\tBlue -\t  Ka: %6.2f Kd: %6.2f Ks: %6.2f\n",
            p->b[0],p->b[1],p->b[2]);

    printf("Red - %6.2f\tGreen - %6.2f\tBlue - %6.2f\n",
            p->red,p->green,p->blue);
}

void print_light(struct light *p) {
    printf("Location -\t %6.2f %6.2f %6.2f\n",
            p->l[0],p->l[1],p->l[2]);

    printf("Brightness -\t r:%6.2f g:%6.2f b:%6.2f\n",
            p->c[0],p->c[1],p->c[2]);
}

void print_symtab() {
    int i;
    for (i=0; i < lastsym;i++) {
        printf("Name: %s\n",symtab[i].name);
        switch (symtab[i].type) {
            case SYM_MATRIX:
                printf("Type: SYM_MATRIX\n");
                print_matrix(symtab[i].s.m);
                break;
            case SYM_CONSTANTS:
                printf("Type: SYM_CONSTANTS\n");
                print_constants(symtab[i].s.c);
                break;
            case SYM_LIGHT:
                printf("Type: SYM_LIGHT\n");
                print_light(symtab[i].s.l);
                break;
            case SYM_VALUE:
                printf("Type: SYM_VALUE\n");
                printf("value: %6.2f\n", symtab[i].s.value);
                break;
            case SYM_STRING:
                printf("Type: SYM_STRING\n");
                printf("Name: %s\n",symtab[i].name);
        }
        printf("\n");
    }
}

SYMTAB *add_symbol(char *name, int type, void *data) {
    SYMTAB *t;

    t = (SYMTAB *)lookup_symbol(name);
    if (t==NULL) {
        if (lastsym >= MAX_SYMBOLS) {
            return NULL;
        }
        t = (SYMTAB *)&(symtab[lastsym]);
        ++lastsym;
        #ifdef DEBUG
        print_debug("Symbol #%d: %s", lastsym, name);
        #endif
    }
    else {
        // If the symbol already exists, return it, and free the argument given.
        // TODO update freeing when implementing other data types
        #ifdef DEBUG
        print_debug("Symbol already exists: %s; "
                    "Freeing add_symbol() argument....", name);
        #endif
        switch (type) {
            case SYM_MATRIX:
                free_matrix((struct matrix *)data);
                break;
            case SYM_CONSTANTS:
                free(data);
                break;
            case SYM_LIGHT:
                free(data);
                break;
            default:
                break;
        }
        return t;
    }

    t->name = (char *)malloc(strlen(name)+1);
    strcpy(t->name,name);
    t->type = type;
    switch (type) {
        case SYM_CONSTANTS:
            t->s.c = (struct constants *)data;
            break;
        case SYM_MATRIX:
            t->s.m = (struct matrix *)data;
            break;
        case SYM_LIGHT:
            t->s.l = (struct light *)data;
            break;
        case SYM_VALUE:
            t->s.value = *(double *)data;
            break;
        case SYM_STRING:
            break;
    }
    return (SYMTAB *)&(symtab[lastsym-1]);
}

SYMTAB *lookup_symbol(char *name) {
    int i;
    for (i=0;i<lastsym;i++) {
        if (!strcmp(name,symtab[i].name)) {
            return (SYMTAB *) &(symtab[i]);
        }
    }
    return (SYMTAB *)NULL;
}

void set_value(SYMTAB *p, double value) {
    p->s.value = value;
}

void free_table() {
    int i;
    #ifdef DEBUG
    print_debug("Freeing symbol table");
    #endif
    for (i = 0; i < lastsym; ++i) {
        // TODO update free_table when implementing other commands
        free(symtab[i].name);
        switch (symtab[i].type) {
            case SYM_MATRIX:
                free_matrix(symtab[i].s.m);
                #ifdef DEBUG
                print_debug("Freeing matrix at %d", i);
                #endif
                break;
            case SYM_CONSTANTS:
                free(symtab[i].s.c);
                #ifdef DEBUG
                print_debug("Freeing constant at %d", i);
                #endif
                break;
            case SYM_LIGHT:
                free(symtab[i].s.l);
                #ifdef DEBUG
                print_debug("Freeing light at %d", i);
                #endif
                break;
            default:
                #ifdef DEBUG
                print_debug("Freeing no additional structures at %d", i);
                #endif
                break;
        }
    }
}
