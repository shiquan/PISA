#include "read_thread.h"
#include "htslib/kstring.h"
#include "read_tags.h"

static char *copy_str(char *a)
{
    if (a == NULL) return NULL;
    return strdup(a);
}
void read_block_clear(struct read_block *b)
{
    int i;
    for (i = 0; i < b->n; ++i) {
        struct read *r = &b->b[i];
        if (r->name) free(r->name);
        if (r->s0) free(r->s0);
        if (r->q0) free(r->q0);
        if (r->s1) free(r->s1);
        if (r->q1) free(r->q1);
    }
    free(b->b);
    if (b->name) free(b->name);
}
struct read_block *read_block_copy(struct read_block *r)
{
    struct read_block *b = malloc(sizeof(*b));
    memset(b, 0, sizeof(*b));
    b->name = copy_str(r->name);
    b->pair_mode = r->pair_mode;
    b->n = r->n;
    b->m = r->n;
    b->b = malloc(sizeof(struct read)*b->n);
    int i;
    for (i = 0; i < b->n; ++i) {
        b->b[i].name = copy_str(r->b[i].name);
        b->b[i].l0 = r->b[i].l0;
        b->b[i].l1 = r->b[i].l1;
        b->b[i].s0 = copy_str(r->b[i].s0);
        b->b[i].q0 = copy_str(r->b[i].q0);
        b->b[i].s1 = copy_str(r->b[i].s1);
        b->b[i].q1 = copy_str(r->b[i].q1);
    }
    return b;
}
void read_block_push(struct read_block *b, struct read *r)
{
    assert(b!=NULL);
    if (b->n == b->m) {
        b->m = b->m == 0 ? 2 : b->m*2;
        b->b = realloc(b->b,b->m*sizeof(struct read));
    }
    b->b[b->n].name = copy_str(r->name);
    b->b[b->n].l0 = r->l0;
    b->b[b->n].l1 = r->l1;
    b->b[b->n].s0 = copy_str(r->s0);
    b->b[b->n].q0 = copy_str(r->q0);
    b->b[b->n].s1 = copy_str(r->s1);
    b->b[b->n].q1 = copy_str(r->q1);
    b->n++;
}
void thread_dat_destroy(struct thread_dat *td)
{
    int i;
    for (i = 0; i < td->n; ++i) {
        struct read_block *rb = &td->rb[i];
        read_block_clear(rb)
    }
    free(td->rb);
    free(td);
}
static char *generate_names(char **names, struct dict *dict)
{
    kstring_t str = {0,0,0};
    int i;
    for (i = 0; i < dict_size(dict); ++i) {
        if (names[i] == NULL) continue;
        kputs(names[i],&str);
    }

    for (i = 0; i < dict_size(dict); ++i) {
        if (names[i] == NULL) continue;
        kstring_t temp = {0,0,0};
        ksprintf(&temp,"|||%s:Z:%s",dict_name(dict,i),names[i]);
        kputs(temp.s, &str);
        free(temp.s);
    }

    return str.s;
}

#define MEM_PER_BLK  100000

static char *name_buf = NULL;
static char *seq_buf = NULL;
static char *qual_buf = NULL;

struct thread_dat *read_thread_dat(FILE *fp, struct dict *tag_dict)
{
    int size = 0;
    char *block_name = NULL;
    char *last_name = NULL;
    int fasta;
    struct thread_dat *td = malloc(sizeof(*td));
    memset(td, 0, sizeof(*td));
    
    if (name_buf != NULL) {
        td->m = 100;
        td->rb = malloc(td->m*sizeof(struct read_block));
        td->n++;
        struct read_block *rb = &td->rb[0];
        memset(rb, 0, sizeof(struct read_block));
        rb->m = 2;
        rb->n = 1;
        rb->b = malloc(sizeof(struct read)*rb->m);
        last_name = name_buf;
        struct read *r = &rb->b[0];
        memset(r, 0, sizeof(struct read));
        r->name  = name_buf;
        r->s0    = seq_buf;
        r->q0    = qual_buf;
        r->l0    = strlen(seq_buf);

        if (tag_dict != NULL) {
            char **tn = fastq_name_pick_tags(name_buf, tag_dict);
            char *bn = generate_names(tn, tag_dict);
            block_name = bn;
            rb->name = block_name;
            if (bn == NULL) error("No tag found at %s", name_buf);
            int i;
            for (i = 0; i < dict_size(tag_dict); ++i) free(tn[i]);
            free(tn);
        }
        // reset buf
        name_buf = NULL;
        seq_buf  = NULL;
        qual_buf = NULL;
    }

    kstring_t name = {0,0,0};
    kstring_t seq  = {0,0,0};
    kstring_t qual = {0,0,0};

    for (;;) {
        name.l = 0;
        seq.l  = 0;
        qual.l = 0;
         
        size++;
        char c = fgetc(fp);        
        if (c == '@') fasta = 0;
        else if (c == '>') fasta = 1;
        else if (c == EOF)  break;
        else error("Unknown input format, only support uncompressed fastq file");

        // kputc(c, &name); // ignore type
        
#define READ_LINE(str) do {\
            c = fgetc(fp); \
            if (c == '\n') break; \
            if (c == '\0') error("Truncated file ?");\
            kputc(c, &str); \
            size++;\
        } while(1) \

        READ_LINE(name);

        READ_LINE(seq);


        if (fasta == 0) {
            c = fgetc(fp);
            if (c != '+') error("Format error. Not fastq file ? '+' vs %c", c);
            for (; c != '\n'; ) c = fgetc(fp); //  qual name line
            READ_LINE(qual);
        }
        char *bn = NULL;
        if (tag_dict != NULL) {
            char **tn = fastq_name_pick_tags(name.s, tag_dict);
            bn = generate_names(tn, tag_dict);
            if (bn == NULL) error("No tag found at %s.", name.s);
        
            int i;
            for (i = 0; i < dict_size(tag_dict); ++i) free(tn[i]);
            free(tn);
        }
        
        struct read_block *b;
        if (block_name == NULL || bn == NULL || strcmp(block_name, bn) != 0) {
            if (size > MEM_PER_BLK) {
                name_buf = strdup(name.s);
                seq_buf = strdup(seq.s);
                qual_buf = qual.s == NULL ? NULL : qual.s;
                free(name.s);
                free(seq.s);
                if (qual.m) free(qual.s);
                return td;
            }
                                        
            block_name = bn;
            if (td->n == td->m) {
                td->m = td->m == 0? 100 : td->m<<1;
                td->rb = realloc(td->rb, sizeof(struct read_block)*td->m);
            }
            b = &td->rb[td->n];
            memset(b, 0, sizeof(struct read_block));
            b->name = bn;
            td->n++;
            last_name = NULL;
        }
        else {
            b = &td->rb[td->n-1];
            free(bn);
        }

        if (last_name == NULL || strcmp(last_name, name.s) != 0) {
            if (b->n == b->m) {
                if (b->m == 0) {
                    b->m = 2;
                    b->b = malloc(2*sizeof(struct read));
                }
                else {
                    b->m = b->m<<1;
                    b->b = realloc(b->b, b->m *sizeof(struct read));
                }
            }
            struct read *read = &b->b[b->n];
            memset(read, 0, sizeof(struct read));
            read->name = strdup(name.s);
            read->s0 = strdup(seq.s);
            read->q0 = qual.s == NULL ? NULL : strdup(qual.s);
            read->l0 = seq.l;
            b->n++;
            last_name = read->name;
        }
        else if (strcmp(last_name, name.s) == 0) {
            struct read *read = &b->b[b->n-1];
            read->s1 = strdup(seq.s);            
            read->q1 = qual.s == NULL ? NULL : strdup(qual.s);
            read->l1 = seq.l;
            last_name = NULL;
            b->pair_mode =1; 
        }
    }

    free(name.s);
    free(seq.s);
    if (qual.m) free(qual.s);

    if (td->n == 0) {
        free(td);
        return NULL;
    }
    return td;
}
