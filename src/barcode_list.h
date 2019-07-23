#ifndef BARCODE_LIST_H
#define BARCODE_LIST_H
#include <stdio.h>
#include <stdlib.h>

struct barcode {
    char *s;
    void *data;
};
struct barcode_list {
    struct barcode *b;
    int n, m;
    void *hash;
};

extern struct barcode_list *barcode_init();
extern int barcode_read(struct barcode_list *lb, const char *fname);
extern int barcode_select(struct barcode_list *lb, char *s);
extern int barcode_push(struct barcode_list *lb, char *s);
extern void barcode_destory(struct barcode_list *b);
#endif
