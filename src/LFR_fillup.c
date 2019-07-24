#include "utils.h"
#include "htslib/khash.h"
#include "htslib/hts.h"
#include "htslib/kstring.h"
#include "htslib/sam.h"
#include "number.h"

struct fill_tags {
    char **tags;
};

KHASH_MAP_INIT_STR(fu, struct fill_tags)

static struct args {
    const char *input_fname;
    const char *output_fname;
    int n_fill;
    char **fill_tags;
    int n_block;
    char **block_tags;

    htsFile *fp;
    htsFile *out;
    bam_hdr_t *hdr;
    int qual_thres;
    int keep_all; // default will filter unannotated reads
} args = {
    .input_fname = NULL,
    .output_fname = NULL,
    .n_fill = 0,
    .fill_tags = NULL,
    .n_block = 0,
    .block_tags = NULL,
    .fp = NULL,
    .out = NULL,
    .hdr = NULL,
    .qual_thres = 20,
    .keep_all = 0,
};

static int parse_args(int argc, char **argv)
{
    int i;
    const char *block_tags = NULL;
    const char *fill_tags = NULL;
    const char *file_thread = NULL;
    for (i = 1; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;
        if (strcmp(a, "-h") == 0) return 1;
        else if (strcmp(a, "-fill") == 0) var = &fill_tags;
        else if (strcmp(a, "-block") == 0 ) var = &block_tags;
        else if (strcmp(a, "-@") == 0) var = &file_thread;
        else if (strcmp(a, "-o") == 0) var = &args.output_fname;
        else if (strcmp(a, "-k") == 0) {
            args.keep_all = 1;
            continue;
        }
        if (var != 0) {
            if (i == argc) error("Miss an argument after %s.", a);
            *var = argv[i++];
            continue;
        }
        if (args.input_fname == NULL) {
            args.input_fname = a;
            continue;
        }
        error("Unknown argument: %s", a);
    }
    
    CHECK_EMPTY(args.output_fname, "-o must be set.");
    CHECK_EMPTY(args.input_fname, "Input bam must be set.");
    CHECK_EMPTY(block_tags, "-block must be set.");
    CHECK_EMPTY(fill_tags, "-fill must be set.");
    kstring_t str = {0,0,0};
    kputs(block_tags, &str);
    int *s = ksplit(&str, ',', &args.n_block);
    assert(args.n_block>0);
    args.block_tags = malloc(sizeof(char*)*args.n_block);
    int j;
    for (j = 0; j < args.n_block; ++j) {
        args.block_tags[j] = strdup(str.s+s[j]);
        if (strlen(args.block_tags[j]) != 2) error("Unknown tag format. Only two character allowed.");
    }
    free(s);
    str.l = 0;
    kputs(fill_tags, &str);
    int *s0 = ksplit(&str, ',', &args.n_fill);
    assert(args.n_fill > 0);
    args.fill_tags = malloc(sizeof(char*)*args.n_fill);
    for (j = 0; j < args.n_fill; ++j) {
        args.fill_tags[j] = strdup(str.s+s0[j]);
        if (strlen(args.fill_tags[j]) != 2) error("Unknown tag format. Only two character allowed.");
    }
    free(s0);
    free(str.s);
    
    int file_th = 5;
    if (file_thread) file_th = str2int((char*)file_thread);
    
    args.fp  = hts_open(args.input_fname, "r");
    CHECK_EMPTY(args.fp, "%s : %s.", args.input_fname, strerror(errno));
    htsFormat type = *hts_get_format(args.fp);
    if (type.format != bam && type.format != sam)
        error("Unsupported input format, only support BAM/SAM/CRAM format.");
    args.hdr = sam_hdr_read(args.fp);
    CHECK_EMPTY(args.hdr, "Failed to open header.");
    //int n_bed = 0;
    args.out = hts_open(args.output_fname, "bw");
    CHECK_EMPTY(args.out, "%s : %s.", args.output_fname, strerror(errno));

    hts_set_threads(args.fp, file_th);
    hts_set_threads(args.out, file_th);
    
    if (sam_hdr_write(args.out, args.hdr)) error("Failed to write SAM header.");
    return 0;
}
static void memory_release()
{
    bam_hdr_destroy(args.hdr);
    sam_close(args.fp);
    sam_close(args.out);
}
struct bam_stack {
    int n, m;
    bam1_t **bam;
    kh_fu_t *dict;
    int n_name, m_name;
    char **names;
};

static struct bam_stack *stack_build()
{
    struct bam_stack *s = malloc(sizeof(*s));
    memset(s, 0, sizeof(*s));
    s->dict = kh_init(fu);
    return s;
}
static void stack_reset(struct bam_stack *s)
{
    int i;
    for (i = 0; i < s->n; ++i)
        if (s->bam[i]) bam_destroy1(s->bam[i]);
    s->n = 0;
    for (i = 0; i < s->n_name; ++i) {
        khint_t k;
        k = kh_get(fu, s->dict, s->names[i]);
        struct fill_tags *t = &kh_val(s->dict, k);
        int j;
        for (j = 0; j < args.n_fill; ++j)
            if (t->tags[j]) free(t->tags[j]);
        free(t->tags);
        kh_del(fu, s->dict,k);
        free(s->names[i]);
    }
    s->n_name = 0;
    // kh_destroy(fu, s->dict);
    // s->dict = kh_init(fu);
}
static void stack_destory(struct bam_stack *s)
{
    int i;
    for (i = 0; i < s->n; ++i)
        if (s->bam[i]) bam_destroy1(s->bam[i]);

    free(s->bam);
    for (i = 0; i < s->n_name; ++i) {
        khint_t k;
        k = kh_get(fu, s->dict, s->names[i]);
        struct fill_tags *t = &kh_val(s->dict, k);
        int j;
        for (j = 0; j < args.n_fill; ++j)
            if (t->tags[j]) free(t->tags[j]);
        free(t->tags);
        // kh_del(fu,s->dict,k);
        free(s->names[i]);
    }
    free(s->names);
    kh_destroy(fu,s->dict);
    free(s);
}
static void stack_finish(struct bam_stack *s)
{
    int i;
    for (i = 0; i < s->n; ++i) {
        // if (s->bam[i] == NULL) continue;
        bam1_t *b = s->bam[i];
        kstring_t str = {0,0,0};
        int j;
        for (j = 0; j < args.n_block; ++j) {
            uint8_t *tag = bam_aux_get(b, args.block_tags[j]);
            if (!tag) {
                if (str.m) free(str.s);
                str.l = 0;
                break;
            }
            kputs((char*)(tag+1), &str);
        }
        if (str.l) {
            khint_t k;
            k = kh_get(fu, s->dict, str.s);
            if (k == kh_end(s->dict)) warnings("%s is not inited.", str.s);
            struct fill_tags *f = &kh_val(s->dict, k);
            int j;
            for (j = 0; j < args.n_fill; ++j) {
                if (f->tags[j] == NULL) continue; // no tag found in this block
                uint8_t *tag = bam_aux_get(b, args.fill_tags[j]);                
                if (!tag) { // empty, fill
                    // debug_print("Append %s to %s", f->tags[j], (char*)b->data);
                    bam_aux_append(b, args.fill_tags[j], 'Z', strlen(f->tags[j])+1, (uint8_t*)f->tags[j]);                    
                }
            }
        }            
    }
}
static void stack_push(struct bam_stack *s, bam1_t *b)
{
    if (s->n == s->m) {
        s->m = s->m == 0 ? 1024 : s->m*2;
        s->bam = realloc(s->bam, s->m*sizeof(void*));
    }
    s->bam[s->n++] = bam_dup1(b);
    kstring_t str = {0,0,0};
    int i;
    for (i = 0; i < args.n_block; ++i) {
        uint8_t *tag = bam_aux_get(b, args.block_tags[i]);
        if (!tag) { // no tag, skip
            if (str.m) free(str.s);
            return;
        }
        kputs((char*)(tag+1), &str);
    }

    khint_t k;
    k = kh_get(fu, s->dict, str.s);
    if (k == kh_end(s->dict)) {
        if (s->n_name == s->m_name) {
            s->m_name += 2;
            s->names = realloc(s->names, s->m_name*sizeof(char*));
        }
        s->names[s->n_name] = strdup(str.s);
        int ret;
        k = kh_put(fu, s->dict, s->names[s->n_name], &ret);
        s->n_name++;
        struct fill_tags *ft = &kh_val(s->dict, k);
        ft->tags = malloc(args.n_fill*sizeof(void*));
        memset(ft->tags,0,args.n_fill*sizeof(void*));
    }
    
    struct fill_tags *ft = &kh_val(s->dict, k);
    for (i = 0; i < args.n_fill; ++i) {
        uint8_t *tag = bam_aux_get(b, args.fill_tags[i]);
        if (!tag) continue;
        char *v = (char*)(tag+1);
        if (ft->tags[i] != NULL) {
            // if (strcmp(ft->tags[i],v) !=0) warnings("Unequal tags in the same block, %s vs %s", ft->tags[i], v);            
        }
        else ft->tags[i] = strdup(v);
    }
    free(str.s);
}
static void stack_write(struct bam_stack *s)
{
    int i;
    for (i = 0; i < s->n; ++i) {
        if (s->bam[i]) {
            if (sam_write1(args.out, args.hdr, s->bam[i]) == -1) error("Failed to write SAM.");
            // bam_destroy1(s->bam[i]);
        }
    }
    stack_reset(s);
}
static int LFR_fillup()
{
    int last_id = -2;
    int ret;
    struct bam_stack *s = stack_build();
    bam1_t *b;
    b = bam_init1();
    while (1) {
        ret = sam_read1(args.fp, args.hdr , b);        
        if (ret < 0) break;
        if (args.keep_all == 0 && b->core.qual < args.qual_thres) continue;
        if (args.keep_all == 0 && b->core.tid < 0) continue;
        if (last_id == -2) last_id = b->core.tid;
        // I could not pre-define how many genes overlapped and how large of each gene, so here only treat each chromsome one by one
        if (b->core.tid != last_id) {
            stack_finish(s);
            stack_write(s);
            last_id = b->core.tid;
        }
        stack_push(s, b);        
    }
    stack_finish(s);
    stack_write(s);
    stack_destory(s);
    bam_destroy1(b);
    return 0;
}
// reads in each block with same BCs will be intepret as fragments from same template
// this program used to fill missed tags for reads from same block.
static int usage()
{
    fprintf(stderr, "LFR_fillup in.bam\n");
    fprintf(stderr, "  -fill         Tags to fill.\n");
    fprintf(stderr, "  -block        Tags to identify each block.\n");
    fprintf(stderr, "  -k            Keep unclassified reads in the output.\n");
    return 1;
}

int LFR_fillup_main(int argc, char **argv)
{
    if (parse_args(argc, argv)) return usage();

    LFR_fillup();

    memory_release();
    
    return 0;
}
