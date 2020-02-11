#include "utils.h"
#include "rld0.h"
#include "rle.h"
#include "mrope.h"
#include "htslib/kstring.h"
#include "htslib/thread_pool.h"
#include "read_thread.h"
#include "read_tags.h"
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
void mag_init_lfr(magopt_t *o)
{
    memset(o, 0, sizeof(magopt_t));
    o->trim_len = 0;
    o->trim_depth = 6;
    
    o->min_elen = 300;
    o->min_ovlp = 0;
    o->min_merge_len = 0;
    o->min_ensr = 1;
    o->min_insr = 0;
    o->min_dratio1 = 0.7;
    
    o->max_bcov = 10.;
    o->max_bfrac = 0.15;
    o->max_bvtx = 64;
    o->max_bdist = 512;
    o->max_bdiff = 50;
}
static void assem_opt_init(fml_opt_t *opt)
{
    opt->n_threads = 1;
    opt->min_asm_ovlp = 20;
    opt->min_merge_len = 0;
    opt->ec_k = -1;
    opt->mag_opt.flag |= MAG_F_AGGRESSIVE;
    opt->mag_opt.flag &= ~MAG_F_POPOPEN;
    opt->min_cnt = 1;
    opt->max_cnt = 2;
    mag_init_lfr(&opt->mag_opt);
};
static struct args {
    const char *input_fname;
    const char *output_fname;
    const char *dist_fname;
    struct dict *tag_dict;
    int pair;
    int n_thread;
    char *last_name;

    // fragment length distribution
    int n_len;
    int *len;
    FILE *fp_dist;
    int mini_overlap;
    fml_opt_t *assem_opt;
    int both_strand;
    uint64_t assem_block;
} args = {
    .input_fname  = NULL,
    .output_fname = NULL,
    .dist_fname   = NULL,
    .tag_dict     = NULL,
    
    .pair         = 0,
    .n_thread     = 1,
    .last_name    = NULL,

    .n_len        = 0,
    .len          = NULL,
    .fp_dist      = NULL,
    
    .mini_overlap = 15,
    .assem_opt    = NULL,
    .both_strand  = 0,
    .assem_block  = 0,
};

extern int assemble_usage();

static int parse_args(int argc, char **argv)
{
    if (argc == 1) return 1;

    int i;
    const char *thread = NULL;
    const char *tags = NULL;
    for (i = 1; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;
        if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return 1;
        else if (strcmp(a, "-t") == 0) var = &thread;
        else if (strcmp(a, "-tag") == 0) var = &tags;
        else if (strcmp(a, "-o") == 0) var = &args.output_fname;
        else if (strcmp(a, "-dis") == 0) var = &args.dist_fname;
        else if (strcmp(a, "-bs") == 0) {
            args.both_strand =1;
            continue;
        }
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
    
    args.tag_dict = str2tag(tags);
    
    if (thread) args.n_thread = str2int((char*)thread);
    if (args.n_thread < 1) args.n_thread = 1;

    args.assem_opt = malloc(sizeof(fml_opt_t));
    assem_opt_init(args.assem_opt);


    if (args.dist_fname) {
        args.fp_dist = fopen(args.dist_fname, "w");
        if (args.fp_dist == NULL) error("%s : %s.", args.dist_fname, strerror(errno));
    }
    return 0;
}
static void memory_release()
{
    dict_destroy(args.tag_dict);
    free(args.assem_opt);
}

struct base_v {
    int l;
    uint8_t *v;
};

static void push_base(struct base_v *v, const uint8_t *s, int l)
{
    uint8_t *e = malloc(l+1);
    int i;
    for (i = 0; i < l; ++i) e[l-i-1] = s[i];
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

static struct rld_t *fmi_gen2(struct base_v *v)
{
    mrope_t *mr;
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

struct read1 {
    kstring_t *name;
    kstring_t *seq;
    kstring_t *qual;
};

struct ret_block {
    int assem_block;
    int n,m;
    struct read1 *a;    
};

static struct ret_block *ret_block_build()
{
    struct ret_block *r = malloc(sizeof(*r));
    memset(r,0,sizeof(*r));
    return r;
}
static kstring_t *kstr_create()
{
    kstring_t *s = malloc(sizeof(*s));
    s->m = s->l = 0;
    s->s = NULL;
    return s;
}
static void kstr_destory(kstring_t *str)
{
    if (str) {
        if (str->m) free(str->s);
        free(str);
    }
}
static void ret_block_destory(struct ret_block *r)
{
    int i;
    for (i = 0; i < r->n; ++i) {
        kstr_destory(r->a[i].name);
        kstr_destory(r->a[i].seq);
        kstr_destory(r->a[i].qual);
    }
    if (r->m) free(r->a);
    free(r);      
}
static void print_utg(fml_utg_t *utg, int n, struct ret_block *r, char *name)
{
    int i;
    for (i = 0; i < n; ++i) {
        const fml_utg_t *u = &utg[i];
        kstring_t *nam  = kstr_create();
        kstring_t *seq  = kstr_create();
        kstring_t *qual = kstr_create();
        
        kputc('@', nam); kputw(i, nam); kputc('_', nam);
        kputs(name,nam);
        kputs("|||SR:i:", nam); kputw(u->nsr, nam);
        kputsn(u->seq, u->len, seq); kputs("", seq);
        kputsn(u->cov, u->len, qual); kputs("", qual);
        if (r->n == r->m) {
            r->m = r->n == 0 ? 2 : r->n*2;
            r->a = realloc(r->a, r->m*sizeof(struct read1));
        }
        struct read1 *a = &r->a[r->n++];
        a->name = nam;
        a->seq = seq;
        a->qual = qual;
    }
}
static void *run_it(void *_d)
{
    struct thread_dat *dat = (struct thread_dat*)_d;
    if (dat->n == 0) {
        thread_dat_destroy(dat);
        return NULL;
    }

    struct ret_block *r = ret_block_build();    
    int i;
    for (i = 0; i < dat->n; ++i) {
        struct read_block *rb = &dat->rb[i];   

        struct base_v *v = rend_bseq(rb);
        
        if (v->l == 0) continue;
        rld_t *e = fmi_gen2(v);
        free(v->v); free(v);
        
        mag_t *g = fml_fmi2mag(args.assem_opt, e);
        mag_g_merge(g, 1, args.mini_overlap);
        mag_g_clean(g, &args.assem_opt->mag_opt);
        mag_g_trim_open(g, &args.assem_opt->mag_opt);

        int n_utg = 0;        
        fml_utg_t *utg = fml_mag2utg(g, &n_utg);
        if (n_utg > 0) {
            print_utg(utg, n_utg, r, rb->name);
            r->assem_block++;
        }
        fml_utg_destroy(n_utg, utg);
    }

    thread_dat_destroy(dat);

    return r;
}

static void write_out(void *s, FILE *out)
{
    if (s == NULL) return;
    struct ret_block *r = (struct ret_block*)s;
    args.assem_block += r->assem_block;
    int i;
    for (i = 0; i < r->n; ++i) {
        struct read1 *a = &r->a[i];
        fputs(a->name->s, out); fputc('\n', out);
        fputs(a->seq->s, out); fputs("\n+\n", out);
        fputs(a->qual->s, out); fputc('\n', out);
        int len = a->seq->l;
        if (len > args.n_len) {
            args.len = args.n_len == 0 ? malloc(len*sizeof(int)) : realloc(args.len, len*sizeof(int));
            for (; args.n_len < len; args.n_len++) args.len[args.n_len] = 0;
        }
        args.len[len-1]++;
    }
    
    ret_block_destory(r);
}

int fastq_assem(int argc, char **argv)
{
    double t_real;
    t_real = realtime();
    
    if (parse_args(argc, argv)) return assemble_usage();

    FILE *fp_in = fopen(args.input_fname, "r");
    if (fp_in == NULL) error("%s : %s.", args.input_fname, strerror(errno));

    FILE *out = stdout;
    if (args.output_fname ) {
        out = fopen(args.output_fname, "w");
        if (out == NULL) error("%s : %s.", args.output_fname, strerror(errno));
    }

    if (args.n_thread == 1) {
        for (;;) {
            struct thread_dat *dat = read_thread_dat(fp_in, args.tag_dict);
            if (dat == NULL) break;
            void *data = run_it(dat);
            write_out(data, out);            
        }
    }
    else {
        hts_tpool *p = hts_tpool_init(args.n_thread);
        hts_tpool_process *q = hts_tpool_process_init(p, args.n_thread*2, 0);
        hts_tpool_result *r;
    
        for (;;) {
            struct thread_dat *dat = read_thread_dat(fp_in, args.tag_dict);
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
    }

    fprintf(stderr, "Assembled block,%"PRIu64"\n", args.assem_block);

    memory_release();

    fclose(fp_in);
    if (args.output_fname)
        fclose(out);
    if (args.fp_dist) {
        int i;
        for (i = 0; i < args.n_len; ++i) {
            if (args.len[i] == 0) continue;
            fprintf(args.fp_dist, "%d,%d\n", i+1, args.len[i]);
        }
        fclose(args.fp_dist);
    }
    LOG_print("Real time: %.3f sec; CPU: %.3f sec", realtime() - t_real, cputime());
    return 0;
}

#ifdef LFR_ASSEM
int main(int argc, char **argv)
{
    return fastq_assem(argc, argv);
}

#endif

