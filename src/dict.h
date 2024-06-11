#ifndef DICT_H
#define DICT_H

#include "utils.h"

struct dict;

struct dict *dict_init();
struct dict *dict_dup(struct dict *D);

void dict_destroy(struct dict *D);

int dict_query(const struct dict *D, char const *key);
int dict_query2(const struct dict *D, char const *key);

int dict_push(struct dict *D, char const *key);
int dict_push1(struct dict *D, char const *key);
int dict_push2(struct dict *D, char const *key, int idx);
int dict_read(struct dict *D, const char *fname, int allow_space);
int dict_read2(struct dict *D, const char *fname, int *val);

char *dict_name(const struct dict *D, int idx);

int dict_size(const struct dict *D);

uint32_t dict_count_sum(const struct dict *D);

uint32_t dict_count(const struct dict *D, int idx);

char **dict_names(struct dict *D);

char *dict_most_likely_key(struct dict *D);

void dict_set_value(struct dict *D);
void *dict_query_value(struct dict *D, int idx);
void *dict_query_value2(struct dict *D, const char *key);
int dict_assign_value(struct dict *D, int idx, void *val);
int dict_delete_value(struct dict *D, int idx);

int dict_del(struct dict *D, const char *key);
int dict_exist(struct dict *D, const char *key);
#endif
