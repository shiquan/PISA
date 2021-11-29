#include "biostring.h"
#include "utils.h"

int *str_split(kstring_t *str, int *_n)
{
    int m=1, n=0;
    int i;
    int *s = malloc(1*sizeof(int));
    s[n++] = 0;
    for (i = 0; i < str->l; ++i) {
        if (str->s[i] == ',' || str->s[i] == ';') {
            if (m == n) {
                m += 2;
                s = realloc(s, sizeof(int)*m);
            }
            s[n++] = i+1;
            str->s[i] = '\0';
        }
    }
    *_n = n;
    return s;
}

