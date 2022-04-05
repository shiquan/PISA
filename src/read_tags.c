// this file contains functions to query and update TAG and values in the name of FASTQ+
// TAG consist of two characters, [A-Za-z][A-Za-z0-9]
// Value is a string or type(on character):string. ACGTN or Z:ACGTN

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
// return pos of tag type, or -1 on not exist
int fname_query_tag(const char *p, const char *tag)
{
    int i, l = strlen(p);
    // optional region of fastq
    for (i = 0; i < l; i++) {
        if (isspace(p[i])) {
            l = i;
            break;
        }
    }
    for (i = 0; i < l-7 && p[i] != '\n' && p[i] != '\0'; ) {
        if (p[i] == '|' && p[i+1] == '|' && p[i+2] == '|') {
            i += 3;
            if (p[i++] != tag[0]) continue;
            if (p[i++] != tag[1]) continue;
            if (p[i++] != ':') continue;
            // i++; // keep the flag
            if (i >= l) break;
            return i;
        }
        i++;        
    }
    return -1;
}

char *fname_pick_tag(const char *p, const char *tag)
{
    int i = fname_query_tag(p, tag);
    if (i < 0) return NULL;
    int l = strlen(p);
    int j = i;
    for (; j<l && p[j] != '|' && !isspace(p[j]);) ++j;
    kstring_t val ={0,0,0};
    kputsn(p+i, j-i, &val);
    kputs("", &val);
    return val.s;
}

char *fname_concat_tag(char *p, const char *tag, const char *val)
{
    kstring_t str = {0,0,0};
    int i, l = strlen(p);
    for (i = 0; i < l && !isspace(p[i]); i++);
    kputsn(p, i, &str); kputs("", &str);
    kputs("|||", &str);
    kputs(tag, &str);
    kputc(':', &str);
    if (val[1] != ':') kputs("Z:", &str); // if no type in the val, treat as 'Z'
    kputs(val, &str);

    if (i < l) kputs(p+i, &str); // other optional fields

    if (str.l > MAX_ID_LENGTH) {
        warnings("Failed to update : length of read name is limited to %d", MAX_ID_LENGTH);
        free(str.s);
        return NULL;
    }
    
    return str.s;
}

char *fname_update_tag(char *p, const char *tag, const char *val)
{
    int i = fname_query_tag(p, tag);    
    if (i < 0) return fname_concat_tag(p, tag, val);
    
    int l = strlen(p);
    int lv = strlen(val);
    int j;
    for (j = i; j < l && !isspace(p[j]); ++j) {
        if (j < l-7 && p[j] == '|' && p[j+1] == '|' && p[j+2] == '|') break;
    }

    if (j - i == lv && strncmp(val, p+i, lv) == 0) return NULL; // already present, no update
    kstring_t str = {0,0,0};
    kputsn(p, i, &str);
    kputc(':', &str);
    kputs(val, &str);
    if (j < l) kputs(p+j, &str);
    if (str.l > MAX_ID_LENGTH) {
        warnings("Failed to update : length of read name is limited to %d", MAX_ID_LENGTH);
        free(str.s);
        return NULL;
    }
    return str.s;
}

char *fname_update_tags(char *p, struct dict *tags, char **vals)
{
    int l = dict_size(tags);
    if (l == 1) return fname_update_tag(p, dict_name(tags, 0), vals[0]);
    
    int i;
    char *new = strdup(p);
    int no_change = 1;
    for (i = 0; i < l; ++i) {
        char *new1 = fname_update_tag(new, dict_name(tags, i), vals[i]);
        if (new1) {
            free(new);
            new = new1;
            no_change = 0;
        }
    }

    if (no_change) {
        free(new);
        return NULL;
    }

    return new;
}

char **fname_pick_tags(const char *p, struct dict *dict)
{
    char **vals = malloc(dict_size(dict)*sizeof(char*));
    memset(vals, 0, dict_size(dict)*sizeof(char*));
    // order value by key
    int i;
    int l = dict_size(dict);
    for (i = 0; i < l; ++i) {
        char *key = dict_name(dict, i);
        char *val = fname_pick_tag(p, key);
        vals[i] = val;
    }
    return vals;
}

char *vals2str(char **vals, int n)
{
    int i;
    kstring_t str = {0,0,0};
    for (i = 0; i < n; ++i) {
        if (i > 0) kputc('_', &str);
        if (vals[i]) kputs(vals[i], &str);
        else {
            if (str.m) free(str.s);
            return NULL;
        }
    }

    for (i = 0; i < str.l; ++i) {
        if (str.s[i] == ':') str.s[i] = '_';
    }
    return str.s;
}

// concate all values into a string
char *fname_tagvalstr(const char *p, struct dict *dict)
{
    char **vals = fname_pick_tags(p, dict);
    int n =  dict_size(dict);
    char *con = vals2str(vals, n);
    int i;
    for (i = 0; i < n; ++i)
        if (vals[i]) free(vals[i]);
    free(vals);

    return con;    
}

