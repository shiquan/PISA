// annotate Gene or Peak to SAM attribution
#include "utils.h"
#include "number.h"
#include "thread_pool.h"
#include "htslib/khash.h"
#include "htslib/kstring.h"
#include "htslib/sam.h"
#include "htslib/khash_str2int.h"
#include "htslib/kseq.h"
#include "htslib/hts.h"
#include "bam_pool.h"
#include "gtf.h"
#include "bed.h"
#include "region_index.h"
#include "dict.h"
#include <zlib.h>

struct read_stat {

    uint64_t reads_in_region;
    uint64_t reads_in_region_diff_strand;
    
    // gtf
    uint64_t reads_in_intergenic;
    uint64_t reads_in_exon;
    uint64_t reads_in_exonintron;
    uint64_t reads_in_intron;
    uint64_t reads_antisense;
    uint64_t reads_ambiguous;
};

static struct args {
    const char *input_fname;
    const char *output_fname;
    const char *bed_fname;
    const char *tag; // attribute in BAM
    const char *gtf_fname;
    const char *report_fname;
    const char *chr_spec_fname;

    const char *group_tag; // export summary information for each group
    
    char **chr_binding;
    const char *btag;
    
    int ignore_strand;
    int splice_consider;
    int intron_consider;
    int n_thread;
    int chunk_size;
    
    htsFile *fp;
    htsFile *out;
    bam_hdr_t *hdr;

    FILE *fp_report;

    // GTF
    struct gtf_spec *G;
    
    // bed
    struct bed_spec *B;
    uint64_t reads_input;
    uint64_t reads_pass_qc;
    int map_qual;

    struct dict *group_stat;

    int debug_mode;
} args = {
    .input_fname     = NULL,
    .output_fname    = NULL,
    .bed_fname       = NULL,    
    .tag             = NULL,    
    .gtf_fname       = NULL,
    .report_fname    = NULL,

    .group_tag       = NULL,
    
    .chr_binding     = NULL,
    .btag            = NULL,
    
    .ignore_strand   = 0,
    .splice_consider = 0,
    .intron_consider = 0,
    .map_qual        = 0,
    .n_thread = 1,
    .chunk_size = 100000,
    .fp              = NULL,
    .out             = NULL,
    .hdr             = NULL,
    .fp_report       = NULL,
    .G               = NULL,    
    .B               = NULL,    
    
    .reads_input     = 0,
    .reads_pass_qc   = 0,
    .group_stat      = 0,
    .debug_mode      = 0
};

static char **chr_binding(const char *fname, bam_hdr_t *hdr)
{
    int n;
    char **list = hts_readlist(fname, 1, &n);
    if (n == 0) error("Empty file.");
    
    char **bind = malloc(n*sizeof(char*));
    struct dict *d = dict_init();
    int i;
    for (i = 0; i < n; ++i) {
        char *s = list[i];
        char *p;
        for (p = s; *p && *p != '\t'; p++);
        if (*p == '\t') *p = '\0';
        int idx = dict_push(d, s);
        if (idx != i) error("%s: Duplicated records ? %s", fname, s);
        s = p +1;
        bind[i] = strdup(s);
        free(list[i]);
    }
    free(list);

    char **ret = malloc(hdr->n_targets*sizeof(char*));
    for (i = 0; i < hdr->n_targets; ++i) {
        int idx = dict_query(d, hdr->target_name[i]);
        if (idx == -1) ret[i] = NULL;
        else {
            ret[i] = bind[idx];
            bind[idx] = NULL;
        }
    }
    
    for (i = 0; i < n; ++i)
        if (bind[i]) free(bind[i]);
    free(bind);
    return ret;
}

// default tags,
// TX for transcript name,
// GN for gene name,
// GX for gene id,
// AN for transcript name but read mapped on antisense,
// RE for region type
//
static char TX_tag[2] = "TX";
static char AN_tag[2] = "AN";
static char GN_tag[2] = "GN";
static char GX_tag[2] = "GX";
static char RE_tag[2] = "RE";

static int parse_args(int argc, char **argv)
{
    int i;
    const char *tags = NULL;
    const char *thread = NULL;
    const char *chunk = NULL;
    const char *file_thread = NULL;
    const char *map_qual = NULL;
    for (i = 1; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;

        // common options
        if (strcmp(a, "-o") == 0 ) var = &args.output_fname;
        else if (strcmp(a, "-report") == 0) var = &args.report_fname;
        else if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return 1;
        else if (strcmp(a, "-t") == 0) var = &thread;
        else if (strcmp(a, "-@") == 0) var = &file_thread;
        else if (strcmp(a, "-chunk") == 0) var = &chunk;
        else if (strcmp(a, "-q") == 0) var = &map_qual;
        else if (strcmp(a, "-debug") == 0) {
            args.debug_mode = 1;
            continue;
        }
        else if (strcmp(a, "-ignore-strand") == 0) {
            args.ignore_strand = 1;
            continue;
        }
        else if (strcmp(a, "-splice-consider") == 0 || strcmp(a, "-splice") == 0) {
            args.splice_consider = 1;
            continue;
        }
        else if (strcmp(a, "-intron") == 0) {
            args.intron_consider = 1;
            continue;
        }
        
        // group options
        else if (strcmp(a, "-group") == 0) var = &args.group_tag;

        // chrom binding
        else if (strcmp(a, "-chr-species") == 0) var = &args.chr_spec_fname;
        else if (strcmp(a, "-btag") == 0) var = &args.btag;
        
        // bed options
        else if (strcmp(a, "-bed") == 0) var = &args.bed_fname;
        else if (strcmp(a, "-tag") == 0) var = &args.tag;

        // gtf options
        else if (strcmp(a, "-gtf") == 0) var = &args.gtf_fname;
        else if (strcmp(a, "-tags") == 0) var = &tags;
        
        if (var != 0) {
            if (i == argc) error("Miss an argument after %s.", a);
            *var = argv[i++];
            continue;
        }
        if (a[0] == '-' && a[1] != '\0') error("Unknown argument, %s",a);
        if (args.input_fname == NULL) {
            args.input_fname = a;
            continue;
        }        
        error("Unknown argument: %s", a);
    }
    
    if (args.bed_fname == NULL && args.gtf_fname == NULL && args.chr_spec_fname == NULL) 
        error("-bed or -gtf or -chr-species must be set.");
    
    CHECK_EMPTY(args.output_fname, "-o must be set.");
    CHECK_EMPTY(args.input_fname, "Input bam must be set.");
    
    if (thread) args.n_thread = str2int((char*)thread);
    if (chunk) args.chunk_size = str2int((char*)chunk);
    if (map_qual) args.map_qual = str2int((char*)map_qual);
    if (args.map_qual < 0) args.map_qual = 0;
    
    int file_th = 1;
    if (file_thread)
        file_th = str2int((char*)file_thread);
    
    args.fp  = hts_open(args.input_fname, "r");
    CHECK_EMPTY(args.fp, "%s : %s.", args.input_fname, strerror(errno));
    htsFormat type = *hts_get_format(args.fp);
    if (type.format != bam && type.format != sam)
        error("Unsupported input format, only support BAM/SAM/CRAM format.");
    args.hdr = sam_hdr_read(args.fp);
    CHECK_EMPTY(args.hdr, "Failed to open header.");
        
    if (args.bed_fname) {
        CHECK_EMPTY(args.tag, "-tag must be set.");
        args.B = bed_read(args.bed_fname);
        if (args.B == 0 || args.B->n == 0) error("Bed is empty.");
    }

    if (args.gtf_fname) {
        LOG_print("Loading GTF.");
        double t_real;
        t_real = realtime();
        
        args.G = gtf_read(args.gtf_fname, 1);
        LOG_print("Load time : %.3f sec", realtime() - t_real);
        
        if (args.G == NULL) error("GTF is empty.");
        if (tags) {
            kstring_t str = {0,0,0};
            kputs(tags, &str);
            if (str.l != 14) error("Bad format of -tags, require five tags and splited by ','.");
            int n;
            int *s = ksplit(&str, ',', &n);
            if (n != 5) error("-tags required five tag names.");
            
            memcpy(TX_tag, str.s+s[0], 2*sizeof(char));
            memcpy(AN_tag, str.s+s[1], 2*sizeof(char));
            memcpy(GN_tag, str.s+s[2], 2*sizeof(char));
            memcpy(GX_tag, str.s+s[3], 2*sizeof(char));
            memcpy(RE_tag, str.s+s[4], 2*sizeof(char));
            free(str.s);
            free(s);
        }
    }

    if (args.chr_spec_fname) {
        if (args.btag == NULL) error("-btag must be set with -chr-species.");
        args.chr_binding = chr_binding(args.chr_spec_fname, args.hdr);
    }

    if (args.report_fname) {
        args.fp_report = fopen(args.report_fname, "w");
        CHECK_EMPTY(args.fp_report, "%s : %s", args.report_fname, strerror(errno));
    }
    else args.fp_report =stderr;

    args.out = hts_open(args.output_fname, "bw");
    CHECK_EMPTY(args.out, "%s : %s.", args.output_fname, strerror(errno));
    if (sam_hdr_write(args.out, args.hdr)) error("Failed to write SAM header.");

    hts_set_threads(args.out, file_th);

    args.group_stat = dict_init();
    int idx;
    idx = dict_push(args.group_stat, "__ALL__");

    dict_set_value(args.group_stat);
    
    struct read_stat *stat = malloc(sizeof(*stat));
    memset(stat, 0, sizeof(*stat));
    dict_assign_value(args.group_stat, idx, stat);
        
    return 0;
}

struct pair {
    int start;
    int end;
};
struct isoform {
    int n;
    struct pair *p;
};
struct isoform *bend_sam_isoform(bam1_t *b)
{
    struct isoform *S = malloc(sizeof(*S));
    memset(S, 0, sizeof(*S));
    int i;
    int start = b->core.pos;
    int l = 0;
    for (i = 0; i < b->core.n_cigar; ++i) {
        int cig = bam_cigar_op(bam_get_cigar(b)[i]);
        int ncig = bam_cigar_oplen(bam_get_cigar(b)[i]);
        if (cig == BAM_CMATCH || cig == BAM_CEQUAL || cig == BAM_CDIFF) {
            l += ncig;
        }
        else if (cig == BAM_CDEL) {
            l += ncig;
        }
        else if (cig == BAM_CREF_SKIP) {
            S->p = realloc(S->p, (S->n+1)*sizeof(struct pair));
            S->p[S->n].start = start +1; // 0 based to 1 based
            S->p[S->n].end = start + l;
            // reset block
            start = start + l + ncig;
            l = 0;
            S->n++;
        }
    } 
    S->p = realloc(S->p, (S->n+1)*sizeof(struct pair));
    S->p[S->n].start = start +1; // 0 based to 1 based
    S->p[S->n].end = start + l;
    S->n++;
    return S;    
}


#include "htslib/thread_pool.h"

struct ret_dat {
    struct bam_pool *p;
    struct dict *group_stat;
    uint64_t reads_input;
    uint64_t reads_pass_qc;
};

static char *exon_type_names[] = {
    "Unknown", "Exon", "Intron", "ExonIntron", "Antisense", "Splice", "Ambiguous", "Intergenic"
};
static char *RE_tags[] = {
    "U", "E", "N", "C", "A", "S", "V", "I",
};
enum exon_type {
    type_unknown = 0,  // unknown type, init state
    type_exon,     // read full covered in exon
    type_intron,        // read full covered in intron, with same strand of gene
    type_exon_intron,   // read cover exon and nearby intron
    type_antisense,     // read map on antisense
    type_splice,        // junction read map two or more exome
    type_ambiguous,     // junction read map to isoform(s) but skip some isoforms between, or map to intron
    type_intergenic,
};

struct trans_type {
    int trans_id;
    enum exon_type type;
};

struct gene_type {
    int gene_id;
    int gene_name;
    enum exon_type type; // main type
    int n, m;
    struct trans_type *a;
};

struct gtf_anno_type {
    enum exon_type type;
    int n, m;
    struct gene_type *a;
};

static void gtf_anno_destroy(struct gtf_anno_type *ann)
{
    int i;
    for (i = 0; i < ann->n; ++i)
        if (ann->a[i].m) free(ann->a[i].a);
    free(ann->a);
    free(ann);
}
static void gtf_anno_print(struct gtf_anno_type *ann, struct gtf_spec const *G)
{
    fprintf(stderr, "Type : %s\n", exon_type_names[ann->type]);
    int i;
    for (i = 0; i < ann->n; ++i) {
        struct gene_type *g = &ann->a[i];
        fprintf(stderr, "Gene : %s, %s\n", dict_name(G->gene_name, g->gene_name),  exon_type_names[g->type]); 
        int j;
        for (j = 0; j < g->n; ++j)
            fprintf(stderr, "  Trans : %s, %s\n", dict_name(G->transcript_id, g->a[j].trans_id), exon_type_names[g->a[j].type]);       
    }
}

// type_exon == type_splice > type_exon_intron > type_ambiguous > type_intron > type_anitisense == type_unknown
static void gene_most_likely_type(struct gene_type *g)
{
    int i;
    for (i = 0; i < g->n; ++i) {
        struct trans_type *t = &g->a[i];
        if (t->type == type_unknown) continue;
        else if (t->type == type_exon) {
            g->type = type_exon;
            return; // exon is already the best hit
        }
        else if (t->type == type_intron) {
            if (g->type == type_unknown) g->type = type_intron; //
        }
        else if (t->type == type_exon_intron) {
            if (g->type == type_unknown) g->type = type_exon_intron; //
            else if (g->type == type_intron) g->type = type_exon_intron; //            
        }
        else if (t->type == type_splice) {
            g->type = type_splice; // same with exon
            return;
        }
        else if (t->type == type_ambiguous) {
            if (g->type == type_unknown) g->type = type_ambiguous; // 
            else if (g->type == type_intron) g->type = type_ambiguous; // should not happen?
        }
    }
    
}
static void gtf_anno_most_likely_type(struct gtf_anno_type *ann)
{
    int i;
    for (i = 0; i < ann->n; ++i) { // refresh type of each gene first
        struct gene_type *g = &ann->a[i];
        gene_most_likely_type(g);
    }
    
    for (i = 0; i < ann->n; ++i) {
        struct gene_type *g = &ann->a[i];
        if (g->type == type_unknown) continue;
        else if (g->type == type_exon) {
            ann->type = type_exon;
            return; // exon is already the best hit
        }
        else if (g->type == type_intron) {
            if (ann->type == type_unknown) ann->type = type_intron; 
        }
        else if (g->type == type_exon_intron) {
            if (ann->type == type_unknown) ann->type = type_exon_intron;
            else if (ann->type == type_intron) ann->type = type_exon_intron; 
        }
        else if (g->type == type_splice) {
            ann->type = type_splice; // same with exon
            return;
        }
        else if (g->type == type_ambiguous) {
            if (ann->type == type_unknown) ann->type = type_ambiguous; // better than unknown
            else if (ann->type == type_intron) ann->type = type_ambiguous;
        }
    }
}

static enum exon_type query_exon(int start, int end, struct gtf const *G, int *exon)
{
    assert(G->type == feature_transcript);
    int j = 0;
    int i;
    for (i = 0; i < G->n_gtf; ++i) {
        // from v0.4, transcript and exon in GTF_Spec struct will be sorted by coordinate
        //struct gtf_lite *g0 = G->strand == 0 ? &G->son[i] : &G->son[G->n_son-i-1];
        struct gtf *g0 = &G->gtf[i];
        if (g0->type != feature_exon) continue;
        j++;
        if (start >= g0->start && end <= g0->end) {
            *exon = j<<2 | (start==g0->start)<<1 | (end == g0->end);            
            return type_exon;
        }

        if (start >= g0->end) continue; // check next exon

        if (end <= g0->start) return type_intron;

        if (start < g0->end && end > g0->end) return type_exon_intron;

        if (start < g0->start && end > g0->start) return type_exon_intron;
    }

    return type_unknown; // out of range
}

// for each transcript, return a type of alignment record
static struct trans_type *gtf_anno_core(struct isoform *S, struct gtf const *g)
{
    struct trans_type *tp = malloc(sizeof(*tp));
    tp->trans_id = g->transcript_id;
    tp->type = type_unknown;
   
    int exon;
    int last_exon = -1;
    int i = 0;
    // linear search
    for (i = 0; i < S->n; ++i) {        
        struct pair *p = &S->p[i];

        enum exon_type t0 = query_exon(p->start, p->end, g, &exon);

        if (t0 == type_unknown) {
            if (tp->type != type_unknown)  tp->type = type_ambiguous; // at least some part of read cover this transcript
            if (i > 0) tp->type = type_ambiguous;
            break; // not covered this exon
        }
        else if (t0 == type_exon) {
            if (tp->type == type_unknown) {
                tp->type = t0;
                last_exon = exon;
                continue;
            }
            else if (tp->type == type_exon || tp->type == type_splice) {
                assert(last_exon != -1);
                if ((exon>>2)-(last_exon >> 2)>1) {  // check the exon number
                    tp->type = type_ambiguous;
                    break;
                }

                if ((last_exon & 0x1) && (exon & 0x2)) { // check the edge
                    tp->type = type_splice;
                    last_exon = exon;
                    continue;
                }
                tp->type = type_ambiguous;
                break;
            }
        }
        else if (t0 == type_intron) {
            if (tp->type == type_unknown) { // intron
                tp->type = t0;
                break;
            }
            else if (tp->type == type_exon || tp->type == type_splice) { // looks like an isoform
                tp->type = type_ambiguous;
                break;                
            }
            else {
                // should not come here
                assert(0);
            }
        }
        else if (t0 == type_exon_intron) {
            tp->type = t0;
            break;
        }
        
    }
    return tp;  
}

// add new trans node to the tree
void gtf_anno_push(struct trans_type *a, struct gtf_anno_type *ann, int gene_id, int gene_name)
{
    if (ann->n == ann->m) {
        ann->m += 5;
        ann->a = realloc(ann->a, ann->m*sizeof(struct gene_type));
    }
    int i;
    for (i = 0; i < ann->n; ++i) {
        struct gene_type *g0 = &ann->a[i];
        if (gene_id == g0->gene_id) {
            if (g0->m == g0->n) {
                g0->m += 5;
                g0->a = realloc(g0->a, sizeof(struct trans_type)*g0->m);
            }
            g0->a[g0->n].trans_id = a->trans_id;
            g0->a[g0->n].type = a->type;
            g0->n++;
            return;
        }
    }

    struct gene_type *g = &ann->a[ann->n];
    ann->n++;
    memset(g, 0, sizeof(*g));
    g->gene_id = gene_id;
    g->gene_name = gene_name;
    g->type = type_unknown;
    g->m = 5;
    g->a = malloc(sizeof(struct trans_type)*g->m);
    g->a[g->n].trans_id = a->trans_id;
    g->a[g->n].type = a->type;
    g->n++;
}
void gtf_anno_string(bam1_t *b, struct gtf_anno_type *ann, struct gtf_spec const *G)
{
    // 
    if (ann->type == type_unknown) return;
    else if (ann->type == type_intron && args.intron_consider == 0) return;  // for default, only annotate intron type 
    else if (ann->type == type_exon_intron && args.splice_consider == 0 && args.intron_consider == 0) return; // 
    else if (ann->type == type_ambiguous) return;

    // only exon or splice come here
    kstring_t gene_name = {0,0,0};
    kstring_t gene_id   = {0,0,0};
    kstring_t trans_id  = {0,0,0};
    int i;
    for (i = 0; i < ann->n; ++i) {
        struct gene_type *g = &ann->a[i];
        if (g->type == ann->type) {
            char *gene = NULL;
            char *id = NULL;
            if (g->gene_name != -1) gene = dict_name(G->gene_name, g->gene_name);
            if (g->gene_id != -1) id = dict_name(G->gene_id, g->gene_id);
            
            if (gene_name.l) {
                kputc(';', &gene_name);
                kputc(';', &gene_id);
                kputc(';', &trans_id);
            }

            if (gene == NULL && id == NULL) error("No gene name or gene id in gtf? %s", (char*)b->data);
            if (gene == NULL) id = gene;
            if (id == NULL) gene = id;
            kputs(gene, &gene_name);
            kputs(id, &gene_id);
            int j;
            int n_trans = 0;
            for (j = 0; j < g->n; ++j) {
                struct trans_type *t = &g->a[j];
                if (t->type == g->type) {                    
                    if (n_trans) kputc(',', &trans_id);
                    char *trans = dict_name(G->transcript_id, t->trans_id);
                    kputs(trans, &trans_id);
                    n_trans++;
                }
            }
        }
    }
    
    if (gene_name.l) {
        bam_aux_append(b, GX_tag, 'Z', gene_id.l+1, (uint8_t*)gene_id.s);
        bam_aux_append(b, GN_tag, 'Z', gene_name.l+1, (uint8_t*)gene_name.s);
        bam_aux_append(b, TX_tag, 'Z', trans_id.l+1, (uint8_t*)trans_id.s);
        free(gene_id.s);
        free(gene_name.s);
        free(trans_id.s);
    }
}

enum exon_type bam_gtf_anno_core(bam1_t *b, struct gtf_spec const *G)
{
    bam_hdr_t *h = args.hdr;
    
    // cleanup all exist tags
    uint8_t *data;
    if ((data = bam_aux_get(b, TX_tag)) != NULL) bam_aux_del(b, data);
    if ((data = bam_aux_get(b, AN_tag)) != NULL) bam_aux_del(b, data);
    if ((data = bam_aux_get(b, GN_tag)) != NULL) bam_aux_del(b, data);
    if ((data = bam_aux_get(b, GX_tag)) != NULL) bam_aux_del(b, data);
    if ((data = bam_aux_get(b, RE_tag)) != NULL) bam_aux_del(b, data);

    bam1_core_t *c;
    c = &b->core;
    
    char *name = h->target_name[c->tid];
    int endpos = bam_endpos(b);

    struct region_itr *itr = gtf_query(G, name, c->pos, endpos);

    // non-overlap, intergenic
    if (itr == NULL) return type_intergenic; // query failed
    if (itr->n == 0) return type_intergenic; // no hit

    enum exon_type ret_type = type_unknown;
    
    struct gtf_anno_type *ann = malloc(sizeof(*ann));
    memset(ann, 0, sizeof(*ann));

    // exon == splice > intron > antisense
    // https://github.com/shiquan/PISA/wiki/4.-Annotate-alignment-records-with-GTF-or-BED

    struct isoform *S = bend_sam_isoform(b);
    
    int antisense = 0;
    int i;
    for (i = 0; i < itr->n; ++i) {
        struct gtf const *g0 = (struct gtf*)itr->rets[i];
        if (g0->start > c->pos+1 || endpos > g0->end) continue; // not fully covered
        
        if (args.ignore_strand == 0) {
            if (b->core.flag & BAM_FREVERSE) {
                if (g0->strand == 0) {
                    antisense = 1;
                    continue;
                }
            }
            else {
                if (g0->strand == 1) {
                    antisense = 1;
                    continue;
                }
            }
        }
        
        int j;
        for (j = 0; j < g0->n_gtf; ++j) {
            struct gtf const *g1 = &g0->gtf[j];
            if (g1->type != feature_transcript) continue;
            struct trans_type *a = gtf_anno_core(S, g1);
            gtf_anno_push(a, ann, g1->gene_id, g1->gene_name);
            free(a);
        }
    }

    // stat type
    gtf_anno_most_likely_type(ann);
    if (ann->type == type_unknown) {
        if (antisense == 1) ret_type = type_antisense;
        else ret_type = type_intergenic; // not fully convered
    }
    else {

        ret_type = ann->type;
                
        gtf_anno_string(b, ann, G);
    }

    if (args.debug_mode) {
        fprintf(stderr, "%s   ", b->data);
        gtf_anno_print(ann, G);
    }
    free(S->p); free(S);
    gtf_anno_destroy(ann);
    region_itr_destroy(itr);

    return ret_type;
}
void bam_gtf_anno(bam1_t *b, struct gtf_spec const *G, struct read_stat *stat)
{
    enum exon_type type = bam_gtf_anno_core(b, G);
    bam_aux_append(b, RE_tag, 'A', 1, (uint8_t*)RE_tags[type]);
    
    if (type == type_exon) stat->reads_in_exon++;
    else if (type == type_splice) stat->reads_in_exon++; // reads cover two exomes
    else if (type == type_intron) stat->reads_in_intron++;
    else if (type == type_exon_intron) stat->reads_in_exonintron++;
    else if (type == type_ambiguous) stat->reads_ambiguous++; // new transcript
    else if (type == type_intergenic) stat->reads_in_intergenic++;
    else if (type == type_antisense) stat->reads_antisense++;
    else error("Unknown type? %s", exon_type_names[type]);
}

void bam_bed_anno(bam1_t *b, struct bed_spec const *B, struct read_stat *stat)
{
    bam_hdr_t *h = args.hdr;
    
    bam1_core_t *c;
    c = &b->core;
    
    // cleanup exist tag
    uint8_t *data;
    if ((data = bam_aux_get(b, args.tag)) != NULL) bam_aux_del(b, data);
    
    char *name = h->target_name[c->tid];
    int endpos = bam_endpos(b);

    struct region_itr *itr = bed_query(args.B, name, c->pos, endpos);
    if (itr == NULL) return; // query failed
    if (itr->n == 0) return; // no hit

    struct dict *val = dict_init();
    int i;
    kstring_t temp = {0,0,0};
    for (i = 0; i < itr->n; ++i) {
        struct bed *bed = (struct bed*)itr->rets[i];
        if (bed->start > endpos || bed->end <= c->pos) continue; // not covered
        temp.l = 0;
        
        if (bed->name == -1) 
            ksprintf(&temp, "%s:%d-%d", dict_name(B->seqname, bed->seqname), bed->start, bed->end);
        else
            kputs(dict_name(B->name, bed->name), &temp);
        
        if (bed->strand == -1)
            stat->reads_in_region++;
        else {
            if (c->flag & BAM_FREVERSE) {
                if (bed->strand == 1) 
                    stat->reads_in_region++;               
                else {
                    stat->reads_in_region_diff_strand++;
                    temp.l = 0; // reset
                }
            }
            else {
                if (bed->strand == 0)
                    stat->reads_in_region++;               
                else {
                    stat->reads_in_region_diff_strand++;
                    temp.l = 0;
                }
            }
        }
        
        if (temp.l)
            dict_push(val, temp.s);
        
    }
    region_itr_destroy(itr);

    if (temp.m) free(temp.s);
    
    if (dict_size(val)) {
        kstring_t str = {0,0,0};
        int i;
        for (i = 0; i < dict_size(val); ++i) {
            if (i) kputc(';', &str);
            kputs(dict_name(val, i), &str);
        }

        bam_aux_append(b, args.tag, 'Z', str.l+1, (uint8_t*)str.s);
        free(str.s);
    }
    dict_destroy(val);
}
void *run_it(void *_d)
{
    bam_hdr_t *h = args.hdr;
    struct ret_dat *dat = malloc(sizeof(struct ret_dat));
    memset(dat, 0, sizeof(*dat));
    dat->p = (struct bam_pool*)_d;
    dat->group_stat = dict_init();

    int idx;
    idx = dict_push(dat->group_stat, "__ALL__");

    dict_set_value(dat->group_stat);
    
    struct read_stat *stat = malloc(sizeof(*stat));
    memset(stat, 0, sizeof(*stat));
    dict_assign_value(dat->group_stat, idx, stat);
    
    int i;
    
    for (i = 0; i < dat->p->n; ++i) {
        bam1_t *b = &dat->p->bam[i];
    
        if (args.group_tag) {
            uint8_t *data;
            if ((data = bam_aux_get(b, args.group_tag)) != NULL) {
                char *tag_val = (char*)(data+1);
                int idx = dict_query(dat->group_stat, tag_val);
                if (idx == -1) {
                    idx = dict_push(dat->group_stat, tag_val);
                    struct read_stat *s = malloc(sizeof(*s));
                    memset(s, 0, sizeof(*s));
                    dict_assign_value(dat->group_stat, idx, s);
                }
                stat = dict_query_value(dat->group_stat, idx);
            }
        }

        dat->reads_input++;
        
        bam1_core_t *c;
        c = &b->core;
        // QC
        if (c->tid <= -1 || c->tid > h->n_targets || (c->flag & BAM_FUNMAP)) continue;
        if (c->qual <= args.map_qual) continue;

        dat->reads_pass_qc++;

        if (args.G) 
            bam_gtf_anno(b, args.G, stat);

        if (args.B)
            bam_bed_anno(b, args.B, stat);
        
        if (args.chr_binding) {
            char *v = args.chr_binding[b->core.tid];
            if (v != NULL)        
                bam_aux_append(b, args.btag, 'Z', strlen(v)+1, (uint8_t*)v);
        }
    }
    return dat;
}

static void write_out(void *_d)
{
    struct ret_dat *dat = (struct ret_dat *)_d;
    int i;
    for (i = 0; i < dat->p->n; ++i) {
        if (sam_write1(args.out, args.hdr, &dat->p->bam[i]) == -1)
            error("Failed to write SAM.");
    }
    
    args.reads_input   += dat->reads_input;
    args.reads_pass_qc += dat->reads_pass_qc;

    for (i = 0; i < dict_size(dat->group_stat); ++i) {
        int idx = dict_query(args.group_stat, dict_name(dat->group_stat, i));
        if (idx == -1) {
            idx = dict_push(args.group_stat, dict_name(dat->group_stat, i));
            struct read_stat *s = malloc(sizeof(*s));
            memset(s, 0, sizeof(*s));
            dict_assign_value(args.group_stat, idx, s);
        }

        struct read_stat *s0 = dict_query_value(args.group_stat, idx);
        struct read_stat *s1 = dict_query_value(dat->group_stat, i);
        s0->reads_in_region += s1->reads_in_region;
        s0->reads_in_region_diff_strand += s1->reads_in_region_diff_strand;
        s0->reads_in_intergenic += s1->reads_in_intergenic;
        s0->reads_in_exon += s1->reads_in_exon;
        s0->reads_in_intron += s1->reads_in_intron;
        s0->reads_antisense += s1->reads_antisense;
        s0->reads_ambiguous += s1->reads_ambiguous;
        s0->reads_in_exonintron += s1->reads_in_exonintron;
    }
    bam_pool_destory(dat->p);

    // free assign memory manually
    for (i = 0; i < dict_size(dat->group_stat); ++i) {
        void *v = dict_query_value(dat->group_stat, i);
        if (v) free(v);
    }

    dict_destroy(dat->group_stat);
    free(dat);
}
void write_report()
{
    if (dict_size(args.group_stat) == 1) {
        struct read_stat *s0 = (struct read_stat*)dict_query_value(args.group_stat, 0);
        fprintf(args.fp_report, "Reads Mapped to Genome (Map Quality > %d),%.1f%%\n", args.map_qual, (float)args.reads_pass_qc/args.reads_input*100);
        
        if (args.B) {
            fprintf(args.fp_report, "Reads Mapped to BED regions / Peaks,%.1f%%\n", (float)s0->reads_in_region/args.reads_pass_qc*100);
            if (s0->reads_in_region_diff_strand)
                fprintf(args.fp_report, "Reads Mapped to BED regions but on diff strand,%.1f%%\n", (float)s0->reads_in_region_diff_strand/args.reads_pass_qc*100);
        }
        if (args.G) {
            fprintf(args.fp_report, "Reads Mapped to Exonic Regions,%.1f%%\n", (float)s0->reads_in_exon/args.reads_pass_qc*100);
            fprintf(args.fp_report, "Reads Mapped to Intronic Regions,%.1f%%\n", (float)s0->reads_in_intron/args.reads_pass_qc*100);
            fprintf(args.fp_report, "Reads Mapped to both Exonic and Intronic Regions,%.1f%%\n", (float)s0->reads_in_exonintron/args.reads_pass_qc*100);
            fprintf(args.fp_report, "Reads Mapped Antisense to Gene,%.1f%%\n", (float)s0->reads_antisense/args.reads_pass_qc*100);
            fprintf(args.fp_report, "Reads Mapped to Intergenic Regions,%.1f%%\n", (float)s0->reads_in_intergenic/args.reads_pass_qc*100);
            fprintf(args.fp_report, "Reads Mapped to Gene but Failed to Interpret Type,%.1f%%\n", (float)s0->reads_ambiguous/args.reads_pass_qc*100);
        }
    }
    else {
    }
}
static void memory_release()
{
    bam_hdr_destroy(args.hdr);
    sam_close(args.fp);
    sam_close(args.out);
    int i;
    for (i = 0; i < dict_size(args.group_stat); ++i) {
        void *v = dict_query_value(args.group_stat, i);
        if (v) free(v);
    }
    dict_destroy(args.group_stat);
    if (args.B) bed_spec_destroy(args.B);
    if (args.G) gtf_destroy(args.G);
    if (args.fp_report != stderr) fclose(args.fp_report);
}

extern int anno_usage();

int bam_anno_attr(int argc, char *argv[])
{
    double t_real;
    t_real = realtime();

    if (parse_args(argc, argv)) return anno_usage();

    if (args.n_thread == 1) {
        for (;;) {
            struct bam_pool *b = bam_pool_create();
            bam_read_pool(b, args.fp, args.hdr, args.chunk_size);
            if (b == NULL) break;
            if (b->n == 0) { free(b->bam); free(b); break; }
            b = run_it(b);
            write_out(b);
        }
    }
    else {
        // multi-thread mode

        hts_tpool *p = hts_tpool_init(args.n_thread);
        hts_tpool_process *q = hts_tpool_process_init(p, args.n_thread*2, 0);
        hts_tpool_result *r;

        for (;;) {
            struct bam_pool *b = bam_pool_create();
            bam_read_pool(b, args.fp, args.hdr, args.chunk_size);
            
            if (b == NULL) break;
            if (b->n == 0) { free(b->bam); free(b); break; }
            
            int block;
            do {
                block = hts_tpool_dispatch2(p, q, run_it, b, 1);
                if ((r = hts_tpool_next_result(q))) {
                    struct bam_pool *d = (struct bam_pool*)hts_tpool_result_data(r);
                    write_out(d);
                    hts_tpool_delete_result(r, 0);
                }
            }
            while (block == -1);
        }

        hts_tpool_process_flush(q);

        while ((r = hts_tpool_next_result(q))) {
            struct bam_pool *d = (struct bam_pool*)hts_tpool_result_data(r);
            write_out(d);
            hts_tpool_delete_result(r, 0);
        }
        hts_tpool_process_destroy(q);
        hts_tpool_destroy(p);
    }

    write_report();
    memory_release();    
    LOG_print("Real time: %.3f sec; CPU: %.3f sec", realtime() - t_real, cputime());
    
    return 0;    
}
