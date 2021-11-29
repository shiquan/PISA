#ifndef BIOSTRING_H
#define  BIOSTRING_H

#include <string.h>
#include "htslib/kstring.h"
int *str_split(kstring_t *str, int *_n);

#endif
