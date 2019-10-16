#ifndef FQ_TAG
#define FQ_TAG
#include <stdio.h>
#include <stdlib.h>
#include "dict.h"

struct dict *str2tag(const char *seq);

extern char **fastq_name_pick_tags(char *name, struct dict *dict);

char *read_name_pick_tag(char *name, const char *tag);

#endif
