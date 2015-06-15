#include "symtab.h"

// TODO use a hash table
SYMTAB symtab[MAX_SYMBOLS];
int lastsym = 0;

void print_constants(struct constants *p) {
    printf("iar: %d\tiag: %d\tiab: %d\n",
           p->iar, p->iag, p->iab);
    printf("kar: %lf\tkag: %lf\tkab: %lf\n",
           p->kar, p->kag, p->kab);
    printf("idr: %d\tidg: %d\tidb: %d\n",
           p->idr, p->idg, p->idb);
    printf("kdr: %lf\tkdg: %lf\tkdb: %lf\n",
           p->kdr, p->kdg, p->kdb);
    printf("isr: %d\tisg: %d\tisb: %d\n",
           p->isr, p->isg, p->isb);
    printf("ksr: %lf\tksg: %lf\tksb: %lf\n",
           p->ksr, p->ksg, p->ksb);
    printf("spec_expt: %d\n", p->spec_expt);
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

    t = lookup_symbol(name);
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
            default:
                #ifdef DEBUG
                print_debug("Freeing no additional structures at %d", i);
                #endif
                break;
        }
    }
}
