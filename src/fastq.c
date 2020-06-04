#include "fastq.h"
#include "number.h"
#include "htslib/kseq.h"
#include "htslib/kstring.h"
#include <zlib.h>

KSEQ_INIT(gzFile, gzread)

struct fastq_pool *fastq_pool_init()
{
    struct fastq_pool *p = malloc(sizeof(*p));
    memset(p, 0, sizeof(*p));
    return p;
}

int trim_read_tail(char *s, int l)
{
    if ( l > 2 && s[l-2] == '/' && (s[l-1] == '1' || s[l-1] == '2')) {
        l -= 2;
        s[l] = '\0';
    }
    return l;
}

void fastq_clean(struct fastq *b)
{
    int i;
    for (i = 0; i < b->n; ++i) {
        free(b->fastq[i].seq);
        if (b->fastq[i].qual) free(b->fastq[i].qual);
    }
    if (b->m_idx) free(b->idx);
    if (b->extend_tags) dict_destroy(b->extand_tags);
    if (b->data) free(b->data);
    if (b->fastq) free(b->fastq);
}
void fastq_destroy(struct fastq *b)
{
    fastq_clean(b);
    free(b);
}
void fastq_pool_destroy(struct fastq_pool *p)
{
    int i;
    for (i = 0; i < p->n; ++i)
        fastq_clean(&p->fastq[i]);
    free(p->fastq);
    free(p);
}

static void fastq_copy(struct fastq *a, struct fastq *b)
{
    // assume a is empty
    memset(a, 0, sizeof(*a));
    
    a->flag = b->flag;
    a->n = b->n;
    a->l_data = b->l_data;
    
    a->extend_tags = dict_copy(b->extend_tags);
    if (a->l_data > 0) {
        a->data = malloc(sizeof(a->l_data));
        memcpy(a->data, b->data, a->l_data);
    }
    
    int i;
    for (i = 0; i < a->n; ++i) {
        a->b[i].length = b->b[i].length;
        a->b[i].seq = strdup(b->b[i].seq);
        a->b[i].qual = b->b[i].qual == NULL ? NULL : strdup(b->b[i].qual);
    }
}
void fastq_pool_push(struct fastq_pool *p, struct fastq *b)
{
    if (p->n == p->m) {
        p->m = p->m == 0 ? 1024 : p->m<<1;
        p->s = realloc(p->s, p->m*sizeof(struct bseq));
    }
    struct fastq *a = &p->s[p->n++];
    fastq_copy(a, b);
}
static char **split_multi_files(const char *fname, int *n)
{
    kstring_t str = {0,0,0};
    kputs(fname, &str);
    int *s = ksplit(&str, ',', n);
    if (*n == 1 || s == 0) {
        free(str.s);
        if (s) free(s);
        return NULL;
    }
    
    char **paths = malloc(*n*sizeof(char*));
    int i;
    for (i = 0; i < *n; ++i) paths[i] = strdup(str.s+s[i]);
    free(str.s);
    free(s);
    return paths;
}

struct ifile {
    char *fname;
    gzFile fp;
    kseq_t *ks;
};

struct input {
    int curr;
    int n;
    struct ifile *in;
};

struct fastq_spec *fastq_spec_init(char **input_fname, int n_file)
{
    struct fastq_spec *h = malloc(sizeof(*h));
    memset(h, 0, sizeof(*h));
    h->n_file = n_file;
    h->input = malloc(h->n_file*sizeof(struct input));
    
    int i;
    int last_n = -1;
    for (i = 0; i < n_file; ++i) { // read 1, 2, 3 or 4
        struct input *in = &h->input[i];
        int n = 0;
        int j;
        char **reads = split_multi_files(input_fname[i], &n);
        in->in = malloc(n*sizeof(struct ifile));
        in->n = n;
        in->curr = 0;
        if (last_n == -1) last_n = n;
        else if (last_n != n) error("Inconsistance fastq input.");
        
        for (j = 0; j < n; ++j) {
            struct ifile *ifile = &in->in[j];
            ifile->fname = reads[j];
            if (n == 1 && n_file == 1) {  // check streaming
                ifile->fp = strcmp(ifile->fname, "-") == 0 ? gzdopen(fileno(stdin), "r"): gzopen(ifile->fname, "r");
            }
            else {
                ifile->fp = gzopen(ifile->fname, "r");
            }
                
            CHECK_EMPTY(ifile->fp, "%s : %s.", ifile->fname, strerror(errno));
            ifile->ks = kseq_init(ifile->fp);
            CHECK_EMPTY(ifile->fp, "Failed to init stream. %s", ifile->fname);
        }
        free(reads);
    }
    
    return h;
}
void fastq_spec_destroy(struct fastq_spec *h)
{
    assert(h->buf == NULL); // make sure buf is empty before destroy 
    int i;
    for (i = 0; i < h->n_file; i++) {
        struct input *in = &h->input[i];
        int j;
        for (j = 0; j < in->n; ++j) {
            struct ifile *ifile = &in->in[j];
            free(ifile->fname);
            if (ifile->fp) gzclose(ifile->fp);
            if (ifile->ks) ks_destroy(ifile->ks);
        }
        free(in->in);
    }
    free(h->input);
    free(h);
}
struct tag_pos {
    int pos
};
void tag_destroy(struct dict *tag)
{
    int i;
    for (i = 0; i < dict_size(tag); ++i) {
        void *data = dict_query_value(tag, i);
        free(data);
    }
    dict_destroy(tag);
}
// this function will change input string
struct dict *fastq_name_parse(char *name, int *len)
{
    int l;
    l = strlen(name);
    if (l == 0) return NULL;

    l = trim_read_tai(name, l);
    
    uint8_t *data = malloc(l);
    struct dict *extend = dict_init();
    dict_set_value(extend);
    int i, j = 0;
    char tag[2];
    for (i = 0; i < l-8; ) {
        if (name[i] == '|' && name[i+1] == '|' && name[i+2] == '|') {
            if (name[i+5] == ':' && name[i+7] == ':') {
                tag[0] = name[i+3];
                tag[1] = name[i+4];
                int idx = dict_push(extend, tag);
                data[j++] = 0;
                struct tag_pos *pos = malloc(sizeof(*pos));
                pos->pos = j;
                dict_assign_value(extend, idx, pos);
                i += 5;
                for (; i < l && (l-i<3 || (name[i] != '|' && name[i+1] != '|' && name[i+2] != '|')); ++i) data[j++] = name[i];
            }
            else {
                i++;
                continue;
            }
        }
        else {
            data[j++] = name[i++];   // read name
        }
    }
    data[j] = 0;
    memset(name, 0, l);
    memcpy(name, data, j);
    *len = j;
    free(data);
    return extend;
}
// add b to a, return
// -1 on inconsistance read name
// 0 on success merged
// 1-N for new tags added
int fastq_check_name(struct fastq *a, char *name)
{
    char *new = strdup(name);
    int len = 0;
    struct dict *extend = fastq_extend(new, &len);
    
    if (a->m_data == 0) {
        a->m_data = len;
        a->data = malloc(a->m_data);
        int i;
        for (i = 0; i < len; ++i) a->data[i] = new[i];
        a->l_data = len;
        free(new);
        a->extend = extend;
        return 0;
    }
    
    if (strcmp(a->data, new) != 0) {
        if (extend) tag_destroy(extend);
        free(new);
        return -1;
    }

    if (a->extend == NULL) {
        a->extend = dict_init();
        dict_set_value(a->extend);
    }
    //
    int new_tag = 0;
    int i;
    for (i = 0; i < dict_size(extend); ++i) {
        int idx = dict_query(a->extend, dict_name(extend, i));
        if (idx == -1) {
            struct tag_pos *new_pos = malloc(sizeof(*new_pos));
            new_pos->pos = a->l_data;
            
            struct tag_pos *pos = dict_query_value(extend, i);
            kstring_t str ={0,0,0};
            kputs(new+pos->pos, &str);
            if (a->l_data+str.l+1 > a->m_data) {
                a->m_data = a->l_data+str.l+1;
                a->data = realloc(a->data, a->m_data);
            }
            int j;
            for (j = 0; j < str.l+1; ++j)
                a->data[++a->l_data] = str.s[j];
            
            dict_set_value(a->extend, idx, new_pos);
            new_tag++;
        }
    }
    return new_tag;
}
struct fastq_pool *fastq_read(struct fastq_spec *fq, int n_record, int max_mem)
{
    assert(fq);
    
    struct fastq_pool *p;
    p = malloc(sizeof(*p));
    memset(p, 0, malloc(sizeof(*p)));
    
    for (;;) { // loop reads
        struct fastq *b = malloc(sizeof(*b));
        memset(b, 0, sizeof(*b));
        b->n = fq->n_file;
        b->b = malloc(b->n*sizeof(struct bseq_core));
        int i;
        for (i = 0; i < fq->n_file; ++i) { // read 1, 2 .., each time read one record
            struct fastq_core *c = &b->b[i];
            memset(c, 0, sizeof(*c));
            
            struct input *in = &fq->input[i];
            struct ifile *file = &in->in[in->curr];
            for (;;) { // if meet the end of file, this loop help reading next file
                int ret;
                ret = kseq_read(file->ks);
                if (ret < 0) {
                    in->curr++; // next handler
                    gzclose(file->fp);
                    kseq_destroy(file->ks);
                    if (in->curr == in->n ) { // already the final one
                        free(b->b);
                        free(b);
                        break;
                    }
                    else {
                        file = &in->in[in->curr];
                        continue;
                    }
                }
                
                ret = fastq_check_name(b, ks->name.s);
                if (ret < 0) error("Inconstance read name. %s vs %s", (char*)b->data, ks->name.s);
                c->length = ks->seq.l;
                c->seq = strdup(ks->seq.s);
                if (ks->qual.l > 0) {
                    if (ks->seq.l != ks->qual.l) error("Inconstance seq and qual length. %s", ks->name.s);
                    c->qual = strdup(ks->qual.s);
                }
                
                break; // always break this loop
            }
        }
        
        if (b == NULL) break; // end of files

        if (fq->input_smart_pair) {
            if (fq->n != 1) error("Only accept 1 input file for smart pairing mode.");
            if (fq->buf) {
                if (strcmp((char*)fp->buf->data, (char*)b->data) == 0) {
                    fastq_merge(b, fp->buf);
                    fastq_destroy(fp->buf);
                }
            }
            else {
                fp->buf = b;
                continue;
            }
        }
        int mem = fastq_pool_push(b, p);
        if (n_record > 0 && p->n >= n_record) break;
        if (max_mem > 0 && max >= mem) break;
    }

    return p;
}
int fastq_pool_dedup(struct fastq_pool *p)
{
}

void *fastq_tag_value(struct fastq *b, const char *tag)
{
    if (b->extend == NULL) return NULL;
    int idx = dict_query(b->extend, tag);
    if (idx == -1) return NULL;
    struct tag_pos *pos = dict_query_value(b->extend, idx);
    return b->data+pos->pos+2;
}
char *fastq_tags(struct fastq *b, struct dict *tags)
{
    int i;
    kstring_t str = {0,0,0};
    for (i = 0; i < dict_size(tags); ++i) {
        char *val = fastq_tag_value(b, dict_name(tags,i));
        if (val) kputs(val, &str);
    }

    return str.s;
}

int fastq_tag_push(struct fastq *b, const char *tag, char type, char *data)
{
    if (b->extend==NULL) {
        b->extend = dict_init();
        dict_set_value(b->extend);
    }

    int idx = dict_query(b->extend, tag);
    if (idx != -1) return -1; // already present

    idx = dict_push(b->extend, tag);
    struct tag_pos *pos = malloc(sizeof(*pos));
    pos->pos = b->l_data;
    kstring_t str = {0,0,0};
    kputc(type, &str);
    kputc(':', &str);
    kputs(data, &str);
    if (b->l_data + str.l +1 < b->m_data) {
        b->m_data = b->l_data + str.l +1;
        b->data = realloc(b->m_data);
    }
    int i;
    for (i = 0; i < str.l+1; ++i) b->data[++b->l_data] = str.s[i];
    return idx;
}

char *fastq_select_seq(struct fastq *b, int rd, int start, int end)
{
    assert(rd < b->n);
    struct fastq_core *c = &b->b[rd];
    if (start > c->length || end > c->length) return NULL; // out of range

    kstring_t str = {0,0,0};
    int len = end -start +1;
    kputsn(c->seq+start-1, len, &str);
    kputs("", &str);
    return str.s;
}
char *fastq_select_qual(struct fastq *b, int rd, int start, int end)
{
    assert(rd < b->n);
    struct fastq_core *c = &b->b[rd];
    if (start > c->length || end > c->length) return NULL; // out of range
    if (c->qual == NULL) return NULL; // for fasta

    kstring_t str = {0,0,0};
    int len = end -start +1;
    kputsn(c->qual+start-1, len, &str);
    kputs("", &str);
    return str.s;
}
int fastq_mean_qual(struct fastq *b)
{
    int qual = 0;
    int length = 0;
    int i, j;
    for (i = 0; i < b->n; ++i) {
        struct fastq_core *c = &b->b[i];
        for (j = 0; j < c->length; ++j) {
            if (c->qual) return -1;
            qual += c->qual[j] - 33;
            length++;
        }
    }
    if (length == 0) return -1;
    return qual/length;
}


/*

void fastq_pool_clean(struct fastq_pool *p)
{
    int i;
    for ( i = 0; i < p->n; ++i ) {
        struct fastq *b = &p->s[i];
        if (b->l0) {
            free(b->n0);
            free(b->s0);
            if (b->q0) free(b->q0);
        }
        if (b->l1) {
            // free(b->n1);
            free(b->s1);
            if (b->q1) free(b->q1);
        }
    }
    if (p->m > 0) free(p->s);
}
void fastq_pool_destroy(struct fastq_pool *p)
{
    fastq_pool_clean(p);
    free(p);
}

static struct fastq_pool *fastq_read_smart(struct fastq_handler *h, int chunk_size)
{
    struct fastq_pool *p = fastq_pool_init();
    int size = 0;
    int ret1= -1;
    do {
        
        ret1 = kseq_read(h->k1);
    
        if (ret1 < 0) { // come to the end of file
            if (h->n_file > 1 && h->curr < h->n_file) {
                gzclose(h->r1);
                kseq_destroy(h->k1);
                h->r1 = gzopen(h->read_1[h->curr], "r");
                h->k1 = kseq_init(h->r1);
                if (kseq_read(h->k1) < 0) error("Empty record ? %s", h->read_1[h->curr]);
            }
            else break;
            h->curr++;
        }

        kseq_t *ks = h->k1;
        
        struct fastq *s;
        if (p->n >= p->m ) {
            p->m = p->m ? p->m*2 : 256;
            p->s = realloc(p->s, p->m*sizeof(struct fastq));
        }
        s = &p->s[p->n];
        memset(s, 0, sizeof(*s));
        trim_read_tail(ks->name.s, ks->name.l);
        s->n0 = strdup(ks->name.s);
        s->s0 = strdup(ks->seq.s);
        s->q0 = ks->qual.l? strdup(ks->qual.s) : 0;
        s->l0 = ks->seq.l;
        size += s->l0;
        if ( kseq_read(ks) < 0 ) error("Truncated input.");

        trim_read_tail(ks->name.s, ks->name.l);
        // s->n1 = strdup(ks->name.s);
        if ( check_name(s->n0, ks->name.s) ) error("Inconsistance paired read names. %s vs %s.", s->n0, ks->name.s);
        s->s1 = strdup(ks->seq.s);
        s->q1 = ks->qual.l? strdup(ks->qual.s) : 0;
        s->l1 = ks->seq.l;
        size += s->l1;
        p->n++;
        if ( size >= chunk_size ) break;
    } while (1);
    if ( p->n == 0 ) {
        fastq_pool_destroy(p);
        return NULL;
   }
    return p;
}
static struct fastq_pool *fastq_read_core(struct fastq_handler *h, int chunk_size, int pe)
{

    // k1 and k2 already load one record when come here

    struct fastq_pool *p = fastq_pool_init();
    int ret1, ret2 = -1;
    
    if ( pe == 0 ) {
        do {
            ret1 = kseq_read(h->k1);
    
            if (ret1 < 0) { // come to the end of file
                if (h->n_file > 1 && h->curr < h->n_file) {
                    gzclose(h->r1);
                    kseq_destroy(h->k1);
                    h->r1 = gzopen(h->read_1[h->curr], "r");
                    h->k1 = kseq_init(h->r1);
                    if (kseq_read(h->k1) < 0) error("Empty record ? %s", h->read_1[h->curr]);
                }
                else break;
                h->curr++;
            }

            if (p->n >= p->m) {
                p->m = p->m ? p->m<<1 : 256;
                p->s = realloc(p->s, p->m*sizeof(struct fastq));
            }
            struct fastq *s = &p->s[p->n];
            kseq_t *k1 = h->k1;
            memset(s, 0, sizeof(*s));
            trim_read_tail(k1->name.s, k1->name.l);
            s->n0 = strdup(k1->name.s);
            s->s0 = strdup(k1->seq.s);
            s->q0 = k1->qual.l? strdup(k1->qual.s) : 0;
            s->l0 = k1->seq.l;
            s->l1 = 0;
            p->n++;
            if ( p->n >= chunk_size ) break;
        }
        while(1);
    }
    else {
        do {            
            ret1 = kseq_read(h->k1);
            ret2 = kseq_read(h->k2);
            if (ret1 < 0) { // come to the end of file
                if (ret2 >=0) error("Inconsistant input fastq records.");
                if (h->n_file > 1 && h->curr < h->n_file) {
                    gzclose(h->r1);
                    kseq_destroy(h->k1);
                    h->r1 = gzopen(h->read_1[h->curr], "r");
                    h->k1 = kseq_init(h->r1);
                    if (kseq_read(h->k1) < 0) error("Empty record ? %s", h->read_1[h->curr]);
                    if (h->r2) {
                        gzclose(h->r2);
                        kseq_destroy(h->k2);
                        h->r2 = gzopen(h->read_2[h->curr], "r");
                        h->k2 = kseq_init(h->r2);
                        if (kseq_read(h->k2) < 0) error("Empty record ? %s", h->read_2[h->curr]);
                    }
                    h->curr++;
                }
                else break;
            }
            kseq_t *k1 = h->k1;
            kseq_t *k2 = h->k2;
            trim_read_tail(k1->name.s, k1->name.l);
            trim_read_tail(k2->name.s, k2->name.l);

            if ( check_name(k1->name.s, k2->name.s) ) error("Inconsistance paired read names. %s vs %s.", k1->name.s, k2->name.s);
            //if ( k1->seq.l != k2->seq.l ) error("Inconsistant PE read length, %s.", k1->name.s);
            
            struct fastq *s;
            if (p->n >= p->m) {
                p->m = p->m ? p->m<<1 : 256;
                p->s = realloc(p->s, p->m*sizeof(struct fastq));
            }
            s = &p->s[p->n];
            memset(s, 0, sizeof(*s));
            s->n0 = strdup(k1->name.s);
            s->s0 = strdup(k1->seq.s);
            s->q0 = k1->qual.l? strdup(k1->qual.s) : 0;
            s->l0 = k1->seq.l;
            //s->n1 = strdup(k2->name.s);
            s->s1 = strdup(k2->seq.s);
            s->q1 = k2->qual.l? strdup(k2->qual.s) : 0;
            s->l1 = k2->seq.l;
            
            p->n++;
            
            if ( p->n >= chunk_size ) break;            
        }
        while(1);
    }
    if ( p->n == 0 ) {
        fastq_pool_destroy(p);
        return NULL;
    }
    return p;
}

struct fastq_handler *fastq_handler_init(const char *r1, const char *r2, int smart, int chunk_size)
{
    struct fastq_handler *h = malloc(sizeof(*h));
    memset(h, 0, sizeof(*h));
    int n1, n2;
    h->read_1 = split_multi_files(r1, &n1);
    if (r2 != NULL) {
        h->read_2 = split_multi_files(r2, &n2);
        if (n1 != n2) error("Unpaired input fastqs.");
    }
    assert(n1 > 0);
    h->n_file = n1;
    h->curr = 1;
    h->smart_pair = smart;
    h->chunk_size = chunk_size;

    if (n1 == 1) {
        h->r1 = strcmp(r1, "-") == 0 ? gzdopen(fileno(stdin), "r") : gzopen(r1, "r");
        if (h->r1 == NULL) error("Failed to open %s : %s.", r1, strerror(errno));
        h->k1 = kseq_init(h->r1);
        if (h->k1 == NULL) error("Failed to init stream. %s", r1);
        
        if (r2) {
            h->r2 = gzopen(r2, "r");    
            if (h->r2 == NULL) error("Failed to open %s: %s.", r2, strerror(errno));
            h->k2 = kseq_init(h->r2);
        }
    }
    else {
        h->r1 = gzopen(h->read_1[0], "r");
        h->k1 = kseq_init(h->r1);
        if (r2) {
            h->r2 = gzopen(h->read_2[0], "r");
            h->k2 = kseq_init(h->r2);
        }
    }
    
    return h;
}
void fastq_handler_destory(struct fastq_handler *h)
{
    kseq_destroy(h->k1);
    gzclose(h->r1);
    if ( h->k2 ) {
        kseq_destroy(h->k2);
        gzclose(h->r2);
    }
    if (h->n_file > 1) {
        int i;
        for (i = 0; i < h->n_file;++i) {
            free(h->read_1[i]);
            if (h->read_2) free(h->read_2[i]);
        }
        free(h->read_1);
        if (h->read_2) free(h->read_2);
    }
    free(h);
}
int fastq_handler_state(struct fastq_handler *h)
{
    if ( h == NULL ) return FH_NOT_ALLOC;
    if ( h->k1 == NULL ) return FH_NOT_INIT;
    if ( h->smart_pair ) return FH_SMART_PAIR;
    if ( h->k2 == NULL ) return FH_SE;
    return FH_PE;    
}
void *fastq_read(void *_h, void *opts)
{
    struct fastq_handler *h = (struct fastq_handler*)_h;

    int state = fastq_handler_state(h);

    struct fastq_pool *b;
    
    switch(state) {
        case FH_SE:
            b = fastq_read_core(h, h->chunk_size, 0);
            break;
            
        case FH_PE:
            b = fastq_read_core(h, h->chunk_size, 1);
            break;
            
        case FH_SMART_PAIR:
            b = fastq_read_smart(h, h->chunk_size);
            break;

        case FH_NOT_ALLOC:
            error("The fastq handler is NOT allocated.");
            break;

        case FH_NOT_INIT:
            error("The fastq handler is NOT inited.");
            break;

        default:
            error("Unknown state");
    }
    
    if (b) b->opts = opts;
    return b;
}

int fastq_process(struct fastq_pool *pool, void *(*func)(void *data))
{
    int i;
    for ( i = 0; i < pool->n; ++i ) {
        struct fastq *s = &pool->s[i];
        s = func(s);
    }
    return 0;
}

#include "htslib/khash.h"

struct hval {
    int idx;
    int qual;
};
struct hvals {
    int n, m;
    struct hval *v;
};
KHASH_MAP_INIT_STR(key, struct hvals*)
//
// Process fastq block...
// stash
void fastq_pool_push(struct fastq *b, struct fastq_pool *p)
{
    assert(p);
    if (p->n == p->m) {
        p->m = p->m == 0 ? 10 : p->m<<1;
        p->s = realloc(p->s, sizeof(struct fastq)*p->m);               
    }
    struct fastq *c = &p->s[p->n++];
    memcpy(c, b, sizeof(struct fastq));
    //debug_print("c %d, b %d", c->l0, b->l0);
    // free(b);
}

// credit to https://github.com/wooorm/levenshtein.c
size_t levenshtein_n(const char *a, const size_t length, const char *b, const size_t bLength) {
  // Shortcut optimizations / degenerate cases.
  if (a == b) {
    return 0;
  }

  if (length == 0) {
    return bLength;
  }

  if (bLength == 0) {
    return length;
  }

  size_t *cache = calloc(length, sizeof(size_t));
  size_t index = 0;
  size_t bIndex = 0;
  size_t distance;
  size_t bDistance;
  size_t result;
  char code;

  // initialize the vector.
  while (index < length) {
    cache[index] = index + 1;
    index++;
  }

  // Loop.
  while (bIndex < bLength) {
    code = b[bIndex];
    result = distance = bIndex++;
    index = SIZE_MAX;

    while (++index < length) {
      bDistance = code == a[index] ? distance : distance + 1;
      distance = cache[index];

      cache[index] = result = distance > result
        ? bDistance > result
          ? result + 1
          : bDistance
        : bDistance > distance
          ? distance + 1
          : bDistance;
    }
  }

  free(cache);

  return result;
}

static char *reverse_seq(char *s, int l)
{
    char *r = malloc(sizeof(char)*l);
    int i;
    for (i = 0; i < l; ++i ) {
        switch(s[i]) {
            case 'A':
                r[l-i-1]='T'; break;
            case 'C':
                r[l-i-1]='G'; break;
            case 'G':
                r[l-i-1]='C'; break;
            case 'T':
                r[l-i-1]='A'; break;
            case 'N':
                r[l-i-1]='N'; break;
            default:
                error("Unknown bases, %s",s);
        }
    }
    return r;
}
static int check_dup(struct fastq *r, struct fastq *q, int strand)
{
    if (strand == -1) error("Unknown strand.");

    kstring_t str = {0,0,0};
    kstring_t str1 = {0,0,0};
    if (r->l1 > 0 && q->l1 > 0) {
        if (strand == 1) {
            int l = q->l1 > r->l0 ? r->l0 : q->l1;
            kputsn(r->s0, l, &str);
            char *rs = reverse_seq(q->s1, q->l1);
            kputsn(rs, l, &str1);
            free(rs);

            // now to read 2

            l = q->l0 > r->l1 ? r->l1 : q->l0;
            kputsn(r->s1+r->l1-l, l, &str);
            char *rs1 = reverse_seq(q->s0, q->l0);
            kputsn(rs1+q->l1-1, l, &str1);
            free(rs1);
        }
        else {
            int l = q->l0 > r->l0 ? r->l0 : q->l0;
            kputsn(q->s0, l, &str);
            kputsn(q->s1, l, &str1);
            l = q->l1 > r->l1 ? r->l1 : q->l1;
            kputsn(r->s1+r->l1-l, l, &str);
            kputsn(q->s1+q->l1-l, l, &str1);
        }
    }
    else {
        if (strand == 1 && q->l0 != r->l0) return 0;
        int l = q->l0 < r->l0 ? q->l0 : r->l0;
        kputsn(r->s0, l, &str);
        kputsn(q->s0, l, &str1);
    }

    kputs("", &str);
    kputs("", &str1);

    assert(str.l == str1.l);
    int score = levenshtein_n(str.s, str.l, str1.s, str1.l);
    if (score > 2) {
        free(str.s);
        free(str1.s);
        return 0;
    }

    free(str.s);
    free(str1.s);
    return 1; // on dup
}
#define DEDUP_SEED 16
// require all sequence should greater than 16bp, and equal length
// the function will destroy all the flags marked before
int fastq_pool_dedup(struct fastq_pool *p)
{
    if (p->n == 1) return 0;
    
    kh_key_t *hash = kh_init(key);
    int n =0, m = 0;
    char **key = NULL;

    int i;
    int j;
    khint_t k;

    kstring_t seed={0,0,0};
    kstring_t rseed = {0,0,0};
    
    for (i = 0; i < p->n; ++i) {
        struct fastq *b = &p->s[i];
        if (b->l0 < DEDUP_SEED) error("Read length is too short. %d", b->l0);
            
        // todo:: convert mark to bits
        b->flag = 0; // reset all the mark
        int qual = 0;
        
        if (b->q0) 
            for (j = 0; j < b->l0; ++j) qual += b->q0[j]-33;
        int strand = -1;
        seed.l = 0; rseed.l = 0; // reset seed string
        kputsn(b->s0, DEDUP_SEED, &seed);
        kputs("", &seed);
        char *rs = b->l1 >0 ? b->s1 : b->s0;
        int l_rs = b->l1 >0 ? b->l1 : b->l0;
        for (j = 0; j < DEDUP_SEED; ++j) {
            switch(rs[l_rs-1-j]) {
                case 'A':
                    kputc('T', &rseed);
                    break;
                case 'C':
                    kputc('G', &rseed);
                    break;
                case 'G':
                    kputc('C', &rseed);
                    break;
                case 'T':
                    kputc('A', &rseed);
                    break;
                case 'N':
                    kputc('N', &rseed);
                    break;
                default:
                    error("Try to reverse unknown base, %s",rs);
            }
            kputs("", &rseed);            
        }
        
        k = kh_get(key, hash, seed.s);
        if (k == kh_end(hash)) {
            // check reverse then
            k = kh_get(key, hash, rseed.s);
            if (k == kh_end(hash)) goto push_to_index; // reverse do not also match
            strand = 1;
            goto check_dup_records;
            
        }
        else {
            strand = 0;
            goto check_dup_records;
        }
        if (0) {
          check_dup_records:
            do {
                struct hvals *v;
                v = kh_val(hash, k);
                for (j = 0; j < v->n; ++j) {
                    struct hval *v1 = &v->v[j];
                    struct fastq *r = &p->s[v1->idx];
                    if (check_dup(r, b, strand)) {
                        //debug_print("%s\t%s\t%d", r->s0, b->s0, v1->idx);
                        if (qual > v1->qual) {
                            // update
                            v1->idx = i;
                            v1->qual = qual;
                            r->flag = FQ_DUP; // set last as dup
                        }
                        else {
                            b->flag = FQ_DUP;
                        }
                        break;
                    }            
                }
                // no found, push to index      
                if (j == v->n) goto push_to_index;
            } while(0);
        }

        if (0) {
          push_to_index:
            // only push forward string to index
            k = kh_get(key, hash, seed.s);
            if (k == kh_end(hash)) {
                if (m == n) {
                    m = m== 0? 10 : m*2;
                    key = realloc(key, m*sizeof(char*));
                }
                key[n] = strdup(seed.s);
                int ret;
                k = kh_put(key, hash, key[n], &ret);
                struct hvals *v = malloc(sizeof(*v));
                memset(v, 0, sizeof(*v));
                v->m = 2; 
                v->v = malloc(sizeof(struct hval)*v->m);
                struct hval *v1 = &v->v[v->n++];
                v1->idx = i;
                v1->qual = qual;
                kh_val(hash, k) = v;                    
                n++; // increase key index
            }
            else {
                struct hvals *v;
                v = kh_val(hash, k);
                if (v->n == v->m) {
                    v->m = v->m*2; // v->m never equal 0
                    v->v = realloc(v->v, v->m*sizeof(struct hval));
                }
                struct hval *v1 = &v->v[v->n++];
                v1->idx = i;
                v1->qual = qual;
            }
        }
    }

    // debug_print("Start debug");
    for (i = 0; i < n; ++i) {
        k = kh_get(key, hash, key[i]);
        if (k == kh_end(hash)) error("Not in the key %s", key[i]);
        struct hvals *v = kh_val(hash, k);
        // debug_print("%s", key[i]);
        free(v->v);
        free(v);
        free(key[i]);
        kh_del(key,hash, k);
    }
    free(seed.s);
    free(rseed.s);
    free(key);
    kh_destroy(key, hash);
    return 0;   
}
/* TODO: improve buffer performance
#define MIN_BUFFER      1     // 1M
#define DEFAULT_BUFFER  1000  // 1G
#define MAX_BUFFER      10000 // 10G

struct fastq_buffer {
    uint32_t n, m;
    struct fastq *s;
    int paired;
    uint32_t l_buf, m_buf;
    uint8_t *buf;
};

struct fastq_handler {
    gzFile read_1;
    gzFile read_2;
    uint8_t *buf;
}
*/
