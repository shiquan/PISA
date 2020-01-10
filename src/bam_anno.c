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
#include "bed_lite.h"
#include "bam_pool.h"
#include "gtf.h"
#include "dict.h"
#include <zlib.h>

static int usage()
{
    fprintf(stderr, "* Annotate bam records with overlapped function regions. Such as gene, trnascript etc.\n");
    fprintf(stderr, "anno_bam -bed peak.bed -tag PK -o anno.bam in.bam\n");
    fprintf(stderr, "anno_bam -gtf genes.gtf -o anno.bam in.bam\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, "  -o               Output bam file.\n");
    fprintf(stderr, "  -q               Mapping quality threshold. [20]\n");
    fprintf(stderr, "  -report          Summary report.\n");
    fprintf(stderr, "  -@               Threads to read and write bam file.\n");
    fprintf(stderr, "\nOptions for BED file :\n");
    fprintf(stderr, "  -bed             Function regions. Three or four columns bed file. Col 4 could be empty or names of this region.\n");
    fprintf(stderr, "  -tag             Attribute tag name. Set with -bed.\n");
    fprintf(stderr, "\nOptions for mixed samples.\n");
    fprintf(stderr, "  -chr-species     Chromosome name and related species binding list.\n");
    fprintf(stderr, "  -btag            Species tag name. Set with -chr-species.\n");
    fprintf(stderr, "\nOptions for GTF file :\n");
    fprintf(stderr, "  -gtf             GTF annotation file. -gtf is conflict with -bed, if set strand will be consider.\n");
    fprintf(stderr, "  -tags            Attribute names. Default is TX,AN,GN,GX,RE.\n");
    fprintf(stderr, "  -ignore-strand   Ignore strand of transcript in GTF. Reads mapped to antisense transcripts will also be count.\n");
    fprintf(stderr, "  -splice-consider Reads covered exon-intron edge will also be count.\n");
    fprintf(stderr, "  -t               Threads to annotate.\n");
    fprintf(stderr, "  -chunk           Chunk size per thread.\n");
    fprintf(stderr, "\nNotice :\n");
    fprintf(stderr, " * For GTF mode, this program will set tags in default, you could also reset them by -tags.\n");
    fprintf(stderr, "   TX : Transcript id.\n");
    fprintf(stderr, "   AN : Same with TX but set only if read mapped to antisense strand of transcript.\n");
    fprintf(stderr, "   GN : Gene name.\n");
    fprintf(stderr, "   GX : Gene ID.\n");
    fprintf(stderr, "   RE : Region type, should E(exon), N(intron)\n");
    return 1;
}

static struct args {
    const char *input_fname;
    const char *output_fname;
    const char *bed_fname;
    const char *tag; // attribute in BAM
    const char *gtf_fname;
    const char *report_fname;
    const char *chr_spec_fname;

    char **chr_binding;
    const char *btag;
    
    int ignore_strand;
    int splice_consider;
    int n_thread;
    int chunk_size;
    
    htsFile *fp;
    htsFile *out;
    bam_hdr_t *hdr;
    FILE *fp_report;
    struct gtf_spec *G;

    uint64_t reads_input;
    uint64_t reads_pass_qc;
    uint64_t reads_in_peak;
    // gtf
    uint64_t reads_in_gene;
    uint64_t reads_in_exon;
    uint64_t reads_in_intron;
    uint64_t reads_antisense;

    int qual_thres;
    // todo: make bedaux more smart
    struct bedaux *B;
    struct bed_chr *last;
    int i_bed;
} args = {
    .input_fname     = NULL,
    .output_fname    = NULL,
    .bed_fname       = NULL,    
    .tag             = NULL,    
    .gtf_fname       = NULL,
    .report_fname    = NULL,

    .chr_binding     = NULL,
    .btag            = NULL,
    .ignore_strand   = 0,
    .splice_consider = 0,
    .n_thread = 4,
    .chunk_size = 100000,
    .fp              = NULL,
    .out             = NULL,
    .hdr             = NULL,
    .fp_report       = NULL,
    .G               = NULL,    
    .B               = NULL,    
    
    .last            = NULL,
    .i_bed           = -1,
    .reads_input     = 0,
    .reads_pass_qc   = 0,
    .reads_in_peak   = 0,
    .reads_in_gene   = 0,
    .reads_in_exon   = 0,
    .reads_in_intron = 0,
    .reads_antisense = 0,

    .qual_thres      = 20,
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
        if (idx != i) error("Duplicated records ? %s", s);        
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

static char TX_tag[2] = "TX";
static char AN_tag[2] = "AN";
static char GN_tag[2] = "GN";
static char GX_tag[2] = "GX";
static char RE_tag[2] = "RE";

static int parse_args(int argc, char **argv)
{
    int i;
    const char *tags = NULL;
    const char *qual = NULL;
    const char *thread = NULL;
    const char *chunk = NULL;
    const char *file_thread = NULL;
    for (i = 1; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;
        if (strcmp(a, "-bed") == 0) var = &args.bed_fname;
        else if (strcmp(a, "-o") == 0 ) var = &args.output_fname;
        else if (strcmp(a, "-report") == 0) var = &args.report_fname;
        else if (strcmp(a, "-tag") == 0) var = &args.tag;
        else if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return 1;
        else if (strcmp(a, "-gtf") == 0) var = &args.gtf_fname;
        else if (strcmp(a, "-tags") == 0) var = &tags;
        else if (strcmp(a, "-q") == 0) var = &qual;
        else if (strcmp(a, "-t") == 0) var = &thread;
        else if (strcmp(a, "-@") == 0) var = &file_thread;
        else if (strcmp(a, "-chunk") == 0) var = &chunk;
        else if (strcmp(a, "-chr-species") == 0) var = &args.chr_spec_fname;
        else if (strcmp(a, "-btag") == 0) var = &args.btag;
        else if (strcmp(a, "-ignore-strand") == 0) {
            args.ignore_strand = 1;
            continue;
        }
        else if (strcmp(a, "-splice-consider") == 0) {
            args.splice_consider = 1;
            continue;
        }
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
    // CHECK_EMPTY(args.bed_fname, "-bed must be set.");

    
    if (args.bed_fname == NULL && args.gtf_fname == NULL && args.chr_spec_fname == NULL) 
        error("-bed or -gtf or -chr-species must be set.");
    
    if (args.bed_fname && args.gtf_fname)
        error("-bed is conflict with -gtf, you can only choose one mode.");
    
    if (args.ignore_strand && args.gtf_fname == NULL)
        error("Only set -ignore-strand with -gtf.");
    
    CHECK_EMPTY(args.output_fname, "-o must be set.");
    CHECK_EMPTY(args.input_fname, "Input bam must be set.");
    if (thread) args.n_thread = str2int((char*)thread);
    if (chunk) args.chunk_size = str2int((char*)chunk);
    if (qual) args.qual_thres = str2int((char*)qual);
    if (args.qual_thres < 0) args.qual_thres = 0; // no filter
    int file_th = 5;
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
    else if (args.gtf_fname) {
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

    hts_set_threads(args.fp, file_th); 
    hts_set_threads(args.out, file_th);

    return 0;
}

int check_is_overlapped_bed(bam_hdr_t *hdr, bam1_t *b, struct bedaux *B)
{
    bam1_core_t *c;
    c = &b->core;
    char *name = hdr->target_name[c->tid];
    if (args.last == NULL || strcmp(B->names[args.last->id], name) != 0) {
        int id = bed_select_chrom(B, name);
        if (id == -1) {
            args.last = NULL;
            return 0;
        }

        args.last = &B->c[id];
        args.i_bed = 0;
    }
    if (args.i_bed == -2) { // out range of bed
        return 0;
    }
    
    for (;;) {
        if (args.i_bed == args.last->n) break;
        if (args.last->b[args.i_bed].end < c->pos+1) args.i_bed++; // iter bed
        else break;
    }

    if (args.i_bed == args.last->n) {
        args.i_bed = -2;
        return 0;
    }
    int end = bam_endpos(b);
    
    if (end < args.last->b[args.i_bed].start) { // read align before region
        return 0;
    }

    kstring_t str = {0,0,0};// name buffer    
    uint8_t *tag = bam_aux_get(b, args.tag);
    if (tag) {
        warnings("%s already present at line %s:%d, skip", args.tag, hdr->target_name[c->tid], c->pos+1);
        return 1;
    }
    int i;
    for (i = args.i_bed; i < args.last->n; ++i) {
        if (end < args.last->b[i].start) break;
        if (str.l) kputc(';', &str);
        if (args.last->b[args.i_bed].name == NULL) {
            ksprintf(&str, "%s:%d-%d", B->names[args.last->id], args.last->b[i].start, args.last->b[i].end);
        }
        else kputs(args.last->b[i].name, &str);
    }
    
    bam_aux_append(b, args.tag, 'Z', str.l+1, (uint8_t*)str.s);
    free(str.s);

    args.reads_in_peak++;
    
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
        else if (cig == BAM_CDEL)            
        {
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

// ret == 0, enclosed in exon
// ret == -1, will reture 3, adjust to -1 thereafter
// ret == 1, (p::start < G::start && p::end > G::start) || (p::start < G::end && p::end > G::end)
// ret == 2, if isoform cover two or more exomes, other cases will not count here. missed isoforms will be count at match_isoform()

static void query_exon(struct pair *p, struct gtf_lite const *G, int *start, int *c, int *ret)
{
    int i;
    int st = -1;
    int ed = -1;

    int ic = *c;
    
    *ret = -2;
    assert(G->type == feature_transcript);
    
    for (i = *start; i < G->n_son; ++i) {
       
        struct gtf_lite const *g0 = G->strand == 0 ? &G->son[i] : &G->son[G->n_son-i-1];
        if (g0->type != feature_exon) continue;
        ic++;
        if (p->start >= g0->start && p->end <= g0->end) {
            *ret = 0;
            *c = ic;
            break;
        }

        if (p->start > g0->end) continue;
        
        if (p->start > g0->start || p->end > g0->start)
            if (st == -1) st = ic; // overlapped exome start

        ed = ic;
        
        if (p->end < g0->start) {
            if (st == -1) {
                *ret = 3;
            }
            break;
        }

        // not break here, because isoform may cover next exome also
        if (p->start < g0->start && p->end >g0->start) *ret = 1;
        if (p->start > g0->end && p->end > g0->end) *ret = 1;
        
    }
    *start = i;

    if (st != -1 && ed != -1 && st != ed) *ret = 1; // cover two or more exomes
    // fprintf(stderr, "%d\t%d\t%d\n",p->start,p->end,*c);
}
//-1 on intron, 0 on match on exon, 1 on overlap exon-intron, 2 on missed isoform at read, 3 on more isoform at read 
int match_isoform(struct isoform *S, struct gtf_lite const *G)
{
    int i;
    int ret = 0;
    int last_c = -1; 
    int c = 0;
    int j = 0;
    for (i = 0; i < S->n; ++i) {

        query_exon(&S->p[i], G, &j, &c, &ret);

        if (S->n == 1 && ret == 3) ret = -1;
        if (last_c == -1) {
            last_c = c;
        }
        else {

            if (c == last_c) { // same exome, should not consider be this transcript. two exomes found at this regions
                ret = 3;
                break;
            }
            /*
            else if (last_c + 1 == c) { // isoform matched, check next then
                continue;
            }
            */
            else if (last_c + 1 < c) { // skip one or more isoform at this transcript
                ret = 2;
                break;
            }
        }
    }
    
    return ret;
}
#include "htslib/thread_pool.h"

struct read_stat {
    uint64_t reads_input;
    uint64_t reads_pass_qc;
    uint64_t reads_in_peak;
    // gtf
    uint64_t reads_in_gene;
    uint64_t reads_in_exon;
    uint64_t reads_in_intron;
    uint64_t reads_antisense;
};

void bam_gtf_anno(struct bam_pool *p, struct read_stat *stat)
{
    struct gtf_itr *itr = gtf_itr_build(args.G);
    
    kstring_t trans = {0,0,0};
    kstring_t genes = {0,0,0};
    kstring_t gene_id = {0,0,0};
    bam_hdr_t *h = args.hdr;
    struct gtf_spec const *G = args.G;
    int i;
    for (i = 0; i < p->n; ++i) {
        stat->reads_input++;
        bam1_t *b = &p->bam[i];
        bam1_core_t *c;
        c = &b->core;
        if (c->tid <= -1 || c->tid > h->n_targets || (c->flag & BAM_FUNMAP)) continue;
        stat->reads_pass_qc++;
        char *name = h->target_name[c->tid];
        int endpos = bam_endpos(b);

        /*
        struct gtf_lite *g;
        int n = 0;
        // todo: improve the algorithm of overlap
        g = gtf_overlap_gene(G, name, c->pos, endpos, &n, 1);
        if (n==0) return 1;
        */
        if ( gtf_query(itr, name, c->pos+1, endpos) != 0 || itr->n == 0) continue;
        
        stat->reads_in_gene++;

        struct isoform *S = bend_sam_isoform(b);
        
        int l;
        int anti = 0;
        int trans_novo = 0;
        int is_intron = 0;
        int is_ambi = 0;
        int splice_overlapped = 0;
        trans.l = genes.l = gene_id.l = 0;
        
        // exon > intron > antisense
        struct gtf_lite const *g = &G->gtf[itr->st];
        
        for (l = 0; l < itr->n; ++l) {            
            struct gtf_lite const *g0 = &g[l];
            if (args.ignore_strand == 0) {
                if (c->flag & BAM_FREVERSE) {
                    if (g0->strand == 0) {
                        anti = 1;
                        continue;
                    }
                }
                else if (g0->strand == 1) {
                    anti = 1;
                    continue;
                }
            }

            if (c->pos+1 < g0->start || endpos > g0->end || c->pos >= g0->end || endpos <= g0->start) {
                is_ambi = 1;
                continue; // not full enclosed in genes
            }
                
            /*
            // for reads overlapped with multiply genes, only count the last one, since the last one is more close with the reads
            g = g + (n-1);
            */
            int i;
            int in_gene = 0;
            // TX
            for (i = 0; i < g->n_son; ++i) {
                struct gtf_lite const *g1 = &g0->son[i];
                if (g1->type != feature_transcript) continue;
                int ret = match_isoform(S, g1);
                if (ret == -2) continue;
                if (ret == 2 || ret == 3) {
                    trans_novo = 1;
                    continue; // not this transcript
                }
                if (ret == -1) {
                    is_intron = 1;
                    continue; // introns will also be filter
                }
                if (ret == 1) { 
                    splice_overlapped = 1;
                    if (args.splice_consider == 0) continue;
                }
                
                trans_novo = 0;
                in_gene = 1;
                if (trans.l) kputc(';', &trans);
                kputs(dict_name(G->transcript_id,g1->transcript_id), &trans);
            }

            if (in_gene) {
                if (genes.l) kputc(';', &genes);
                kputs(dict_name(G->gene_name,g0->gene_name), &genes);
                if (gene_id.l) kputc(';', &gene_id);
                kputs(dict_name(G->gene_id,g0->gene_id), &gene_id);
                
                anti = 0; // reset antisense flag
            } 
        }
        // GN_tag
        // char *gene = G->gene_name->name[g->gene_name];
        // l = strlen(gene);
        if (trans.l) { // match case
            bam_aux_append(b, TX_tag, 'Z', trans.l+1, (uint8_t*)trans.s);
            bam_aux_append(b, RE_tag, 'A', 1, (uint8_t*)"E");
            bam_aux_append(b, GN_tag, 'Z', genes.l+1, (uint8_t*)genes.s);
            bam_aux_append(b, GX_tag, 'Z', gene_id.l+1, (uint8_t*)gene_id.s);
            stat->reads_in_exon++;
        }
        else if (anti) { // antisense
            bam_aux_append(b, RE_tag, 'A', 1, (uint8_t*)"A");
            stat->reads_antisense++;
        }
        else if (trans_novo) {
            kputs(dict_name(G->gene_name,g->gene_name), &genes);
            bam_aux_append(b, GN_tag, 'Z', genes.l+1, (uint8_t*)genes.s);
            bam_aux_append(b, TX_tag, 'Z', 8, (uint8_t*)"UNKNOWN");            
        }
        else if (is_intron) {
            bam_aux_append(b, RE_tag, 'A', 1, (uint8_t*)"I");
            stat->reads_in_intron++;
        }
        else if (is_ambi) {
            bam_aux_append(b, RE_tag, 'A', 1, (uint8_t*)"U");
        }
        else if (splice_overlapped) {
            bam_aux_append(b, RE_tag, 'A', 1, (uint8_t*)"S");
        }
        
        free(S->p);
        free(S);
    }
    if (trans.m) free(trans.s);
    if (genes.m) free(genes.s);
    if (gene_id.m) free(gene_id.s);

    gtf_itr_destory(itr);
}

struct ret_dat {
    struct bam_pool *p;
    struct read_stat *stat;
};
void *run_it(void *_d)
{    
    struct ret_dat *dat = malloc(sizeof(struct ret_dat));
    dat->p = (struct bam_pool*)_d;
    dat->stat = malloc(sizeof(struct read_stat));
    memset(dat->stat, 0, sizeof(struct read_stat));
    
    if (args.G) bam_gtf_anno(dat->p, dat->stat);

    if (args.chr_binding) {
        int i;
        for (i = 0; i < dat->p->n; ++i) {
            bam1_t *b = &dat->p->bam[i];
            //debug_print("%d",b->core.tid);
            if (b->core.tid == -1) continue;
            char *v = args.chr_binding[b->core.tid];
            if (v == NULL) continue;
            //debug_print("%s",v);
            bam_aux_append(b, args.btag, 'Z', strlen(v)+1, (uint8_t*)v);
        }
    }
    
    return dat;
}

static void write_out(void *_d)
{
    struct ret_dat *dat = (struct ret_dat *)_d;
    int i;
    for (i = 0; i < dat->p->n; ++i) 
        if (sam_write1(args.out, args.hdr, &dat->p->bam[i]) == -1) error("Failed to write SAM.");

    args.reads_in_peak += dat->stat->reads_in_peak;
    args.reads_input += dat->stat->reads_input;
    args.reads_pass_qc += dat->stat->reads_pass_qc;
    args.reads_in_gene += dat->stat->reads_in_gene;
    args.reads_in_exon += dat->stat->reads_in_exon;
    args.reads_in_intron += dat->stat->reads_in_intron;
    args.reads_antisense += dat->stat->reads_antisense;

    bam_pool_destory(dat->p);
    free(dat->stat);
    free(dat);
}
void write_report()
{
    if (args.B) {
        fprintf(args.fp_report, "Reads Mapped Confidently to Genome,%.1f%%\n", (float)args.reads_pass_qc/args.reads_input*100);
        fprintf(args.fp_report, "Reads Mapped Confidently to Peaks,%.1f%%\n", (float)args.reads_in_peak/args.reads_pass_qc*100);
    }
    if (args.G) {
        fprintf(args.fp_report, "Reads Mapped Confidently to Genome,%.1f%%\n", (float)args.reads_pass_qc/args.reads_input*100);
        fprintf(args.fp_report, "Reads Mapped Confidently to Gene,%.1f%%\n", (float)args.reads_in_gene/args.reads_input*100);
        fprintf(args.fp_report, "Reads Mapped Confidently to Exonic Regions,%.1f%%\n", (float)args.reads_in_exon/args.reads_input*100);
        fprintf(args.fp_report, "Reads Mapped Confidently to Intronic Regions,%.1f%%\n", (float)args.reads_in_intron/args.reads_input*100);
        fprintf(args.fp_report, "Reads Mapped Antisense to Gene,%.1f%%\n", (float)args.reads_antisense/args.reads_input*100);
    }
}
static void memory_release()
{
    bam_hdr_destroy(args.hdr);
    sam_close(args.fp);
    sam_close(args.out);
    if (args.B) bed_destroy(args.B);
    else if (args.G) gtf_destory(args.G);
    if (args.fp_report != stderr) fclose(args.fp_report);
}

int bam_anno_attr(int argc, char *argv[])
{
    double t_real;
    t_real = realtime();

    if (parse_args(argc, argv)) return usage();

    if (args.B) { // todo: support multi-threads
        bam1_t *b;
        int ret;
        b = bam_init1();
        while ((ret = sam_read1(args.fp, args.hdr, b)) >= 0) {
            args.reads_input++;
            // todo: QC?
            if (b->core.qual < args.qual_thres) continue;
            args.reads_pass_qc++;
            check_is_overlapped_bed(args.hdr, b, args.B); 
            if (sam_write1(args.out, args.hdr, b) == -1) error("Failed to write SAM.");
        }
        bam_destroy1(b);
        if (ret != -1) warnings("Truncated file?");   
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
