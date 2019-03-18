#ifndef BARCODE_LIST_H
#define BARCODE_LIST_H
#include <stdio.h>
#include <stdlib.h>

struct barcode {
    char *s;
    void *data;
};
struct lbarcode {
    struct barcode *b;
    int n, m;
    void *hash;
};

extern struct lbarcode *barcode_init();
extern int barcode_read(struct lbarcode *lb, const char *fname);
extern int barcode_select(struct lbarcode *lb, char *s);
extern int barcode_push(struct lbarcode *lb, char *s);
extern void barcode_destory(struct lbarcode *b);
#endif
