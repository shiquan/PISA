#include "read_tags.h"
#include "htslib/kstring.h"
#include <ctype.h>

struct dict *str2tag(const char *seq)
{
    kstring_t str = {0,0,0};
    kputs(seq, &str);
    int n;
    int *s = ksplit(&str, ',', &n);
    struct dict *dict = dict_init();
    int i;
    for (i = 0; i <n; ++i) dict_push(dict,str.s+s[i]);
    free(s); free(str.s);
    return dict;
}
char *read_name_pick_tag(char *p, const char *tag)
{
    int i, l;
    l = strlen(p);
    for (i = 0; i < l-7 && p[i] != '\n' && p[i] != '\0'; ) {
        if (p[i] == '|' && p[i+1] == '|' && p[i+2] == '|') {
            i += 3;
            if (p[i++] != tag[0]) continue;
            if (p[i++] != tag[1]) continue;
            if (p[i++] != ':') continue;
            i++; // skip flag
            if (i >= l) break;
            if (p[i++] != ':') continue;
            int j = i;
            for (;i<l && p[i] != '|' && !isspace(p[i]);) ++i;
            kstring_t val ={0,0,0};
            kputsn(p+j, i-j, &val);
            kputs("", &val);
            return val.s;
        }
        else {
            i++;
        }
    }
    return NULL;
}
char **fastq_name_pick_tags(char *p, struct dict *dict)
{
    int l;
    l = strlen(p);
    char key[3];
    char **vals = malloc(dict_size(dict)*sizeof(char*));
    memset(vals, 0, dict_size(dict)*sizeof(char*));
    int i;
    for (i = 0; i < l-7; ) {
        if (p[i] == '|' && p[i+1] == '|' && p[i+2] == '|') {
            i += 3;
            key[0] = p[i++];
            key[1] = p[i++];
            key[2] = '\0';
            if (p[i++] != ':') continue;

            int id = dict_query(dict, key);
            if (id == -1) continue;
            // todo: check type
            i++; // skip flag
            if (i >= l) break;
            if (p[i++] != ':') continue;
            int j = i;
            for (;i<l && p[i] != '|';) ++i;
            kstring_t val ={0,0,0};
            kputsn(p+j, i-j, &val);
            kputs("", &val);
            vals[id] = val.s;
        }
        else {
            i++;
        }
    }

    return vals;    
}

