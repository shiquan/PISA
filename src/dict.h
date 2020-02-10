#ifndef DICT_H
#define DICT_H

#include "utils.h"

struct dict;

struct dict *dict_init();

void dict_destroy(struct dict *D);

int dict_query(const struct dict *D, char *key);

int dict_push(struct dict *D, char *key);

int dict_read(struct dict *D, const char *fname);

char *dict_name(const struct dict *D, int idx);

int dict_size(const struct dict *D);

uint32_t dict_count_sum(const struct dict *D);

uint32_t dict_count(const struct dict *D, int idx);

char **dict_names(struct dict *D);

char *dict_most_likely_key(struct dict *D);
#endif
