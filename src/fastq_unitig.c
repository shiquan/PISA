#include "utils.h"
#include "rld0.h"
#include "rle.h"
#include "mrope.h"
#include "htslib/kstring.h"
#include "htslib/thread_pool.h"
#include "fastq.h"
#include "ksw.h"
#include "fml.h"
#include "number.h"
#include "mag.h"

extern int ksa_sa(const unsigned char *T, int *SA, int n, int k);
extern int ksa_bwt(unsigned char *T, int n, int k);

static unsigned char nt6_tab[256] = {
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5
};

static void revcomp6(uint8_t *s, int l)
{
	int i;
	for (i = 0; i < l>>1; ++i) {
		int tmp = s[l-1-i];
		tmp = (tmp >= 1 && tmp <= 4)? 5 - tmp : tmp;
		s[l-1-i] = (s[i] >= 1 && s[i] <= 4)? 5 - s[i] : s[i];
		s[i] = tmp;
	}
	if (l&1) s[i] = (s[i] >= 1 && s[i] <= 4)? 5 - s[i] : s[i];
}
static uint8_t *enc_str(char *s, int l)
{
    uint8_t *e = malloc(l+1);
    int i;
    for (i = 0; i < l; ++i) e[i] = nt6_tab[(int)s[i]]; 
    e[l] = 0;
    return e;
}
static void assem_opt_init(fml_opt_t *opt)
{
    fml_opt_init(opt);
    opt->min_asm_ovlp = 10;
    opt->ec_k = -1;
    opt->mag_opt.flag = MAG_F_NO_SIMPL | MAG_F_AGGRESSIVE;
};
static struct args {
    const char *input_fname;
    const char *output_fname;
    struct dict *tag_dict;
    int pair;
    int min_ovlp;
    int n_thread;
    char *last_name;

    fml_opt_t *assem_opt;
    int l_seed;
    uint8_t *seed;
    uint64_t assem_block;
} args = {
    .input_fname = NULL,
    .output_fname = NULL,
    .tag_dict = NULL,

    .pair = 0,
    .min_ovlp = 10,
    .n_thread = 1,
    .last_name = NULL,

    .assem_opt = NULL,
    .l_seed = 0,
    .seed = NULL,
    .assem_block = 0,

};

static int usage()
{
    fprintf(stderr, "fastq_unitig in.fq\n");
    fprintf(stderr, "   -t         Threads.\n");
    fprintf(stderr, "   -o         Output fastq.\n");
    fprintf(stderr, "   -tag       Tags of read block.\n");
    fprintf(stderr, "   -ss        Seed sequence for scLFR library. Usually be fixed adaptor sequences.\n");
    fprintf(stderr, "   -p         Input fastq is smart paired.\n");
    return 1;
}

static int parse_args(int argc, char **argv)
{
    if (argc == 1) return 1;

    int i;
    const char *thread = NULL;
    const char *tags = NULL;
    const char *seed = NULL;
    for (i = 1; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;
        if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return 1;
        else if (strcmp(a, "-t") == 0) var = &thread;
        else if (strcmp(a, "-tag") == 0) var = &tags;
        else if (strcmp(a, "-o") == 0) var = &args.output_fname;
        else if (strcmp(a, "-ss") == 0) var = &seed;
        else if (strcmp(a, "-p") == 0) {
            args.pair = 1;
            continue;
        }

        if (var != 0) {
            *var = argv[i++];
            continue;
        }

        if (args.input_fname == NULL) {
            args.input_fname = a;
            continue;
        }

        error("Unknown argument : %s", a);
    }
    if (args.input_fname == NULL ) error("No input fastq specified.");
    if (tags == NULL) error("-tag must be set.");
    if (seed) { //error("Seed sequence for scLFR must be set with -ss.");
        args.l_seed = strlen(seed);
        if (args.l_seed < 6) error("Seed is too short, need >= 6bp.");
        args.seed = enc_str((char*)seed, args.l_seed);
    }
    
    kstring_t str = {0,0,0};
    kputs(tags, &str);
    int n;
    int *s = ksplit(&str, ',', &n);
    args.tag_dict = dict_init();

    for (i = 0; i <n; ++i) dict_push(args.tag_dict,str.s+s[i]);
    free(s); free(str.s);

    if (thread) args.n_thread = str2int((char*)thread);
    if (args.n_thread < 1) args.n_thread = 1;

    args.assem_opt = malloc(sizeof(fml_opt_t));
    assem_opt_init(args.assem_opt);
    return 0;
}
static void memory_release()
{
    dict_destroy(args.tag_dict);
}
struct read {
    char *name;
    int l0, l1;
    char *s0;
    char *q0;
    char *s1;
    char *q1;
};

struct read_block {
    char *name;
    struct read *b;
    int n, m;
    int pair_mode;
};

struct thread_dat {
    int n, m;
    struct read_block *rb;
};

void thread_dat_destroy(struct thread_dat *td)
{
    int i;
    for (i = 0; i < td->n; ++i) {
        struct read_block *rb = &td->rb[i];
        free(rb->name);
        int j;
        for (j = 0; j < rb->n; ++j) {
            struct read *read = &rb->b[j];
            free(read->name);
            free(read->s0);
            if (read->q0) free(read->q0);
            if (read->s1) free(read->s1);
            if (read->q1) free(read->q1);
        }
        free(rb->b);
    }
    free(td->rb);
    free(td);
}
static char *generate_names(char **names)
{
    kstring_t str = {0,0,0};
    int i;
    for (i = 0; i < dict_size(args.tag_dict); ++i) {
        if (names[i] == NULL) continue;
        kputs(names[i],&str);
    }

    for (i = 0; i < dict_size(args.tag_dict); ++i) {
        if (names[i] == NULL) continue;
        kstring_t temp = {0,0,0};
        ksprintf(&temp,"|||%s:Z:%s",dict_name(args.tag_dict,i),names[i]);
        kputs(temp.s, &str);
        free(temp.s);
    }

    return str.s;
}

#define MEM_PER_BLK  512000000

static char *name_buf = NULL;
static char *seq_buf = NULL;
static char *qual_buf = NULL;

struct thread_dat *read_thread_dat(FILE *fp)
{
    int size = 0;
    char *block_name = NULL;
    char *last_name = NULL;
    int fasta;
    struct thread_dat *td = malloc(sizeof(*td));
    memset(td, 0, sizeof(*td));
    
    if (name_buf != NULL) {
        td->m = 10000;
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
        char **tn = fastq_name_pick_tags(name_buf, args.tag_dict);
        char *bn = generate_names(tn);

        block_name = bn;
        rb->name = block_name;
        if (bn == NULL) error("No tag found at %s", name_buf);
        int i;
        for (i = 0; i < dict_size(args.tag_dict); ++i) free(tn[i]);
        free(tn);
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

        kputc(c, &name);
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

        char **tn = fastq_name_pick_tags(name.s, args.tag_dict);
        char *bn = generate_names(tn);
        if (bn == NULL) error("No tag found at %s.", name.s);
        
        int i;
        for (i = 0; i < dict_size(args.tag_dict); ++i) free(tn[i]);
        free(tn);

        struct read_block *b;
        if (block_name == NULL || strcmp(block_name, bn) != 0) {
            if (size > MEM_PER_BLK) {
                name_buf = strdup(name.s);
                seq_buf = strdup(seq.s);
                qual_buf = qual.s == NULL ? NULL : qual.s;
                return td;
            }
                                        
            block_name = bn;
            if (td->n == td->m) {
                td->m = td->m == 0? 10000 : td->m<<1;
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

struct base_v {
    int l;
    uint8_t *v;
};

static void push_base(struct base_v *v, const uint8_t *s, int l)
{
    uint8_t *e = malloc(l+1);
    // memcpy(e, s, l+1);
    int i;
    for (i = 0; i < l; ++i) e[l-i-1] = s[i];
    // for (i = 0; i < l; ++i) e[i] = s[i];
    e[l] = 0;
    v->v = realloc(v->v, v->l+(l+1)*2);    
    memcpy(v->v+v->l, e, l+1);
    v->l += l+1;
    revcomp6(e, l);    
    memcpy(v->v + v->l, e, l+1);
    v->l += l+1;
    free(e);
}

static void push_str_base(struct base_v *v, char *s)
{
    if (s == NULL) return;
    int l;
    l = strlen(s);
    int i;
    for (i = 0; i < l; ++i)
        if (nt6_tab[(int)s[i]] == 5) return;
    uint8_t *e = enc_str(s, l);
    push_base(v, e, l);
    free(e);
}

static struct base_v *rend_bseq(struct read_block *b)
{
    struct base_v *v = malloc(sizeof(*v));
    v->l = 0; v->v = NULL;

    int i;
    for (i = 0; i < b->n; i++) {        
        push_str_base(v, b->b[i].s0);
        if (b->b[i].l1)
            push_str_base(v, b->b[i].s1);
    }
    
    return v;
}

static int bwt_gen(int asize, int64_t l, uint8_t *s)
{
    if (l <= INT32_MAX) return ksa_bwt(s, l, asize);
    else error("Do not support so many reads.");
}

static rld_t *bwt_enc(int asize, int sbits, int64_t l, const uint8_t *s)
{
	int c;
	int64_t i, k;
	rlditr_t itr;
	rld_t *e;

	e = rld_init(asize, sbits);
	rld_itr_init(e, &itr, 0);
	k = 1; c = s[0];
	for (i = 1; i < l; ++i) {
		if (s[i] != c) {
			rld_enc(e, &itr, k, c);
			c = s[i];
			k = 1;
		} else ++k;
	}
	rld_enc(e, &itr, k, c);
	rld_enc_finish(e, &itr);
	return e;
}
/*
static int check_polyTs(const rld_t *e, int len)
{
    uint64_t k,l, ok, ol;
    int i;
    k = e->cnt[4];
    l = e->cnt[5]-1;
    for (i = len -2; i>=0; --i) {
        rld_rank21(e, k, l, 4, &ok, &ol);
        k = e->cnt[4]+ok;
        l = e->cnt[4]+ol;
        if (k >= l) break;
    }
    if ( k >= l ) return 1;
    return 0;
}
*/

int bwt_search_id(const rld_t *e, const char *s, int len)
{
    debug_print("s: %s", s);
    uint64_t ok, ol, k, l;
    int c = nt6_tab[(int)s[len-1]];
    k = e->cnt[c];
    l = e->cnt[c+1]-1;
    
    int i;
    // for (i = 1; i < len ; ++i) {
    for (i = len-2; i >= 10; --i) {
        c = nt6_tab[(int)s[i]];
        if (c == 5) return -1;
        // c = c== 0 ? 0 : 5-c;
        rld_rank21(e, k, l, c, &ok, &ol);
        //debug_print("%c\t%lld\t%lld","$ACGTN"[c],ok,ol);
        k = e->cnt[c] + ok;
        l = e->cnt[c] + ol;
        if (k==l) break;        
    }

    // fprintf(stderr, "\n");
    assert(k<=l);
    // debug_print("K: %d, L: %d", k, l);
    if (k != l-1) return -1;
    k = l;
    uint64_t k6[6];
    while (1) {
        int c = rld_rank1a(e, k, k6);
        fputc("$ACGTN"[c],stderr);
        if (c == 0)  fputc('\n', stderr);
        k = e->cnt[c] + k6[c];
        //debug_print("l: %lld, cnt: %d, occ: %d, %c", k, e->cnt[c], k6[c], "$ACGTN"[c]);;
        if (c == 0) return k;
    }

    return -1;
}
static int bwt_backward_search(const rld_t *e, int len, const uint8_t *str, uint64_t *sa_beg, uint64_t *sa_end)
{
	uint64_t k, l, ok, ol;
	int i, c;
        c = str[len-1];
        
	k = e->cnt[c];
        l = e->cnt[c + 1]-1;
        for (i = len-2; i >= 0; --i) {
            c = str[i];
            rld_rank21(e, k, l, c, &ok, &ol);
            k = e->cnt[c] + ok;
            l = e->cnt[c] + ol;
            // debug_print("%d\t%c\t%d\t%d",i,"$ACGTN"[c], k, l);
            if (k == l) break;
	}
        // if (i < len) return -1;
	// if (k == l) return 0;
	*sa_beg = k; *sa_end = l;
	return l - k;
}

static rld_t *bwt_build_core(rld_t *e0, int asize, int sbits, int64_t l, uint8_t *s)
{
	rld_t *e;
        bwt_gen(asize, l, s);
        // putchar('\n');                
        e = bwt_enc(asize, sbits, l, s);
	return e;
}

rld_t *bwt_build(struct base_v *v)
{
    return bwt_build_core(0, 6, 3, v->l, v->v);
}
int64_t bwt_retrieve(const rld_t *e, uint64_t x, kstring_t *s)
{
    uint64_t k = x, *ok;
    ok = alloca(8 * e->asize);
    s->l = 0;
    while (1) {
        int c = rld_rank1a(e, k, ok);
        k = e->cnt[c] + ok[c] -1;
        if (c == 0) return k;
        kputc(c, s);
    }
}

struct rld_t *fmi_gen2(struct base_v *v)
{
    mrope_t *mr;
    //kstring_t str = {0,0,0};
    mritr_t itr;
    rlditr_t di;
    const uint8_t *block;
    rld_t *e = 0;
    mr = mr_init(ROPE_DEF_MAX_NODES, ROPE_DEF_BLOCK_LEN, MR_SO_RCLO);
    mr_insert_multi(mr, v->l, v->v, 0);
    e = rld_init(6, 3);
    rld_itr_init(e, &di, 0);
    mr_itr_first(mr, &itr, 1);
    while ((block = mr_itr_next_block(&itr)) != 0) {
        const uint8_t *q = block + 2, *end = block + 2 + *rle_nptr(block);
        while (q < end) {
            int c = 0;
            int64_t l;
            rle_dec1(q, c, l);
            rld_enc(e, &di, l, c);
        }
    }
    rld_enc_finish(e, &di);
    
    mr_destroy(mr);
    return e;
}

struct ret_block {
    int assem_block;
    char *s;
};

struct ret_block *ret_block_build()
{
    struct ret_block *r = malloc(sizeof(*r));
    memset(r, 0, sizeof(*r));
    return r;
}
// After construct unitigs from reads, we map reads back to the unitigs. By using
// the paired reads, we connect unitigs into a scaffold
struct kmer_idx {
    struct dict *dict;
    int m, n;
    struct idx {
        int n, m;
        struct offset {
            int strand;
            int idx;
            int pos;
        } *offset;
    } *idx;
    int n_z, m_z;
    char **z; // cache original sequences
};

struct kmer_idx *kmer_idx_build()
{
    struct kmer_idx *i = malloc(sizeof(*i));
    memset(i, 0, sizeof(*i));
    i->dict = dict_init();
    return i;
}

void kmer_idx_destroy(struct kmer_idx *idx)
{
    int i;
    for (i = 0; i < idx->n; i++)
        if (idx->idx[i].m)
            free(idx->idx[i].offset);
    free(idx->idx);
    dict_destroy(idx->dict);
    for (i = 0; i < idx->n_z; ++i) free(idx->z[i]);
    free(idx->z);
    free(idx);
}
#define SEED_LEN 30

void kmer_idx_add(struct kmer_idx *idx, char *s, int l)
{
    if (l < SEED_LEN) error("Read is too short. Require at least %d bases.", SEED_LEN);

    char seed[SEED_LEN];
    int strand;
    for (strand = 0; strand < 2; ++strand) {
        char *z = strdup(s);
        if (strand == 1) {
            int i;
            for (i = 0; i < l>>1; ++i) {
                char c = "$ACGTN"[5-nt6_tab[(int)z[i]]];
                z[i] = "$ACGTN"[5-nt6_tab[(int)z[l-i-1]]];
                z[l-i-1] = c;
            }
            if (l&1) {
                i = l/2+1;
                z[i] = "$ACGTN"[5-nt6_tab[(int)z[i]]];
            }
        }
        if (idx->n_z == idx->m_z) {
            idx->m_z = idx->m_z == 0 ? 10 : idx->m_z<<1;
            idx->z = realloc(idx->z, idx->m_z*sizeof(char*));
        }
        idx->z[idx->n_z] = z;
        int i;
        for (i = 0; i < l - SEED_LEN; ++i) {        
            memcpy(seed, z+i, SEED_LEN);

            int id = dict_query(idx->dict, seed);
            if (id == -1) {
                id = dict_push(idx->dict, seed);
                if (id >= idx->m) {
                    idx->m = id+1;
                    idx->idx = realloc(idx->idx, idx->m*sizeof(struct idx));
                    for ( ; idx->n < idx->m; idx->n++) 
                        memset(&idx->idx[idx->n], 0, sizeof(struct idx));                    
                }
            }
            struct idx *d = &idx->idx[id];
            if (d->n == d->m) {
                d->m = d->m == 0 ? 1 : d->m<<1;
                d->offset = realloc(d->offset, sizeof(struct offset)*d->m);
                d->offset[d->n].strand = strand;
                d->offset[d->n].idx = idx->n_z;
                d->offset[d->n].pos = i;
                d->n++;
            }
        }
        idx->n_z++;
    }
}

int kmer_query_id(struct kmer_idx *idx, char *s, int l)
{
    if (l < SEED_LEN) error("Read is too short. Require at least %d bases.", SEED_LEN);
    char seed[SEED_LEN];
    memcpy(seed, s, SEED_LEN);
    int off = dict_query(idx->dict, seed);
    if (off == -1) return -1;
    int ret = -1;
    int i;
    for (i = 0; i < idx->idx[off].n; ++i) {
        struct offset *offset = &idx->idx[off].offset[i];
        char *z = idx->z[offset->idx];
        s = s + SEED_LEN;
        z = z + offset->pos + SEED_LEN;
        int k;
        for (k = 0; k < l-SEED_LEN;  ++k)
            if (s[k] != z[k]) break;

        if (k < l - SEED_LEN) continue;
        if (ret != -1) return -1; // too much hits
        ret = offset->idx/2;
    }

    //debug_print("ret : %d",ret);
    return ret;
}
static char *remap_reads_scaf(struct read_block *rb, mag_t *g)
{
    struct kmer_idx *idx = kmer_idx_build();
    int i;
    kstring_t str = {0,0,0};
    for (i = 0; i < g->v.n; ++i) {
        magv_t *gv = &g->v.a[i];
        if (gv->len < 0) continue;
        str.l = 0;
        int j;
        for (j = 0; j < gv->len; ++j) kputc("$ACGTN"[(uint8_t)gv->seq[j]], &str);
        kmer_idx_add(idx, str.s, str.l);        
    }
    int *bidx = malloc(g->v.n*sizeof(int));
    memset(bidx, 0, g->v.n*sizeof(int));
    int max_bidx = 1;
    
    for (i = 0; i < rb->n; ++i) {
        if (rb->b[i].s1 == NULL) continue;
        int id1 = kmer_query_id(idx, rb->b[i].s0, rb->b[i].l0);
        int id2 = kmer_query_id(idx, rb->b[i].s1, rb->b[i].l1);
        if (id1 == -1 || id2 == -1) continue;
        
        if (bidx[id1] != 0 && bidx[id2] != 0) {
            if (bidx[id1] != bidx[id2]) {
                int c = bidx[id2];
                int j;
                for (j =0; j < g->v.n; ++j)
                    if (bidx[j] == c) bidx[j] = bidx[id1];
            }
        }
        else {
            if (bidx[id1] != 0) bidx[id2] = bidx[id1];
            else if (bidx[id2] != 0) bidx[id1] = bidx[id2];
            else bidx[id1] = bidx[id2] = max_bidx++;
        }
    }
    kmer_idx_destroy(idx);
    str.l = 0;

    for (i = 0; i < g->v.n; ++i) {
        magv_t *v = &g->v.a[i];
        if (v->len < 0) continue;
        
        ksprintf(&str, "@%d_%s|||SR:i:%d",i, rb->name,v->nsr);

        if (bidx[i] != 0) ksprintf(&str, "|||RB:Z:%d", bidx[i]);
        kputc('\n', &str);
        int j;
        for (j = 0; j < v->len; ++j) kputc("$ACGTN"[(int)v->seq[j]], &str);
        kputs("\n+\n", &str);
        kputsn(v->cov, v->len, &str);
        kputc('\n', &str);            
    }
    free(bidx);
    return str.s;
}

static char *remap_reads(struct read_block *rb, mag_t *g)
{
    kstring_t s = {0,0,0};
    struct base_v *v = malloc(sizeof(*v));
    memset(v, 0, sizeof(*v));
    
    int i;
    for (i = 0; i < g->v.n; ++i) {
        magv_t *v = &g->v.a[i];
        if (v->len < 0) continue;
        
        ksprintf(&s, "@%d_%s|||SR:i:%d\n",i, rb->name,v->nsr);
        int j;
        for (j = 0; j < v->len; ++j) kputc("$ACGTN"[(int)v->seq[j]], &s);
        kputs("\n+\n", &s);
        kputsn(v->cov, v->len, &s);
        kputc('\n', &s);            
    }
    return s.s;
}

static void *run_it(void *_d)
{
    struct thread_dat *dat = (struct thread_dat*)_d;
    if (dat->n == 0) {
        thread_dat_destroy(dat);
        return NULL;
    }

    struct ret_block *r = ret_block_build();
    kstring_t str = {0,0,0};
    // r->all_block = dat->n;
    int i;
    for (i = 0; i < dat->n; ++i) {
        struct read_block *rb = &dat->rb[i];   

        struct base_v *v = rend_bseq(rb);
        
        if (v->l == 0) continue;
        
        rld_t *e = fmi_gen2(v);
        free(v->v); free(v);
        
        if (args.l_seed) {        
            int n;
            uint64_t st = 0, ed = 0;
            n = bwt_backward_search(e, args.l_seed, args.seed, &st, &ed);
            if (n == 0) continue;
        }

        mag_t *g = fml_fmi2mag(args.assem_opt, e);

        char *s = NULL;
        if (rb->pair_mode)
            s = remap_reads_scaf(rb, g);
        else
            s = remap_reads(rb, g);

        if (s) {
            kputs(s, &str);
            free(s);
            r->assem_block++;
        }
        
        mag_g_destroy(g);
    }

    thread_dat_destroy(dat);
    
    if (str.l == 0) {
        free(r);
        return NULL;
    }
    r->s = str.s;
    return r;
}

static void write_out(void *s, FILE *out)
{
    if (s == NULL) return;
    struct ret_block *r = (struct ret_block*)s;
    fputs(r->s, out);
    args.assem_block += r->assem_block;
    free(r->s);
    free(r);
}

int fastq_unitig(int argc, char **argv)
{
    double t_real;
    t_real = realtime();
    
    if (parse_args(argc, argv)) return usage();

    FILE *fp_in = fopen(args.input_fname, "r");
    if (fp_in == NULL) error("%s : %s.", args.input_fname, strerror(errno));

    
    FILE *out = stdout;
    if (args.output_fname ) {
        out = fopen(args.output_fname, "w");
        if (out == NULL) error("%s : %s.", args.output_fname, strerror(errno));
    }
    
    hts_tpool *p = hts_tpool_init(args.n_thread);
    hts_tpool_process *q = hts_tpool_process_init(p, args.n_thread*2, 0);
    hts_tpool_result *r;
    
    for (;;) {
        struct thread_dat *dat = read_thread_dat(fp_in);
        if (dat == NULL) break;

        int block;
        do {
            block = hts_tpool_dispatch2(p, q, run_it, dat, 1);
            if ((r = hts_tpool_next_result(q))) {
                void *s = hts_tpool_result_data(r);
                write_out(s, out);
                hts_tpool_delete_result(r, 0);
            }
        }
        while (block == -1);
    }

    hts_tpool_process_flush(q);
    
    while ((r = hts_tpool_next_result(q))) {
        void *s = hts_tpool_result_data(r);
        write_out(s, out);
        hts_tpool_delete_result(r, 0);
    }
    
    hts_tpool_process_destroy(q);
    hts_tpool_destroy(p);
    fclose(out);
    fprintf(stderr, "Assembled block,%"PRIu64"\n", args.assem_block);

    memory_release();

    fclose(fp_in);
    fclose(out);
    LOG_print("Real time: %.3f sec; CPU: %.3f sec", realtime() - t_real, cputime());
    return 0;
}

/*
int main(int argc, char **argv)
{
    return assem(argc, argv);
}
*/
