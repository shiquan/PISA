#ifndef DICT_H
#define DICT_H

#include "utils.h"

struct dict;

struct dict *dict_init();

void dict_destroy(struct dict *D);

int dict_query(const struct dict *D, char const *key);
int dict_queryInt(const struct dict *D, int key);

int dict_push(struct dict *D, char const *key);

int dict_read(struct dict *D, const char *fname);

char *dict_name(const struct dict *D, int idx);
int dict_nameInt(const struct dict *D, int idx);

int dict_size(const struct dict *D);

uint32_t dict_count_sum(const struct dict *D);

uint32_t dict_count(const struct dict *D, int idx);

char **dict_names(struct dict *D);

char *dict_most_likely_key(struct dict *D);

void dict_set_value(struct dict *D);
void *dict_query_value(struct dict *D, int idx);
void *dict_query_value2(struct dict *D, const char *key);
void *dict_query_valueInt(struct dict *D, int key);
int dict_assign_value(struct dict *D, int idx, void *val);

int dict_pushInt(struct dict *D, int key);
int dict_pushInt1(struct dict *D, int key);


#endif
