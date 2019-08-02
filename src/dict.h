#ifndef DICT_H
#define DICT_H

#include "utils.h"

struct dict;

struct dict *dict_init();

void dict_destroy(struct dict *D);

int dict_query(struct dict *D, char *key);

int dict_push(struct dict *D, char *key);

int dict_read(struct dict *D, const char *fname);

char *dict_name(struct dict *D, int idx);

int dict_size(struct dict *D);

uint32_t dict_count_sum(struct dict *D);

uint32_t dict_count(struct dict *D, int idx);

#endif
