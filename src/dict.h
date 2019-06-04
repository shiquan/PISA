#ifndef DICT_H
#define DICT_H

struct dict {
    int n, m;
    char **name;
    void *dict;
};

struct dict *dict_init();
void dict_destory(struct dict *d);
int dict_query(struct dict *d, char *name);
int dict_push(struct dict *d, char *name);

#endif
