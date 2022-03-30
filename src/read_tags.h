#ifndef FQ_TAG
#define FQ_TAG
#include <stdio.h>
#include <stdlib.h>
#include "dict.h"
#include "fastq.h"

struct dict *str2tag(const char *seq);
int fname_query_tag(const char *p, const char *tag);

char *fname_pick_tag(const char *p, const char *tag);
char *fname_concat_tag(char *p, const char *tag, const char *val);
char *fname_update_tag(char *p, const char *tag, const char *val);
char *fname_update_tags(char *p, struct dict *tags, char **vals);
char **fname_pick_tags(const char *p, struct dict *dict);
char *fname_tagvalstr(const char *p, struct dict *dict);

char *vals2str(char **vals, int n);

#endif
