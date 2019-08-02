// Count all reads and reads with predefined tag for each cell barcode in the bam file
#include "utils.h"
#include "number.h"
#include "htslib/khash.h"
#include "htslib/kstring.h"
#include "htslib/sam.h"
#include "htslib/kseq.h"
#include "dict.h"
#include <zlib.h>

static int usage()
{
    fprintf(stderr, "BamCount in.bam\n");
    fprintf(stderr, "Options :\n");
    fprintf(stderr, "    -cb          Cell Barcode, or other tag used for each individual.\n");
    fprintf(stderr, "    -list        Cell barcode white list.\n");
    fprintf(stderr, "    -tags        Tags to count.\n");
    fprintf(stderr, "    -dedup       Deduplicate the atrributes in each tag.\n");
    fprintf(stderr, "    -group       Group tag, count all tags for each group seperately.\n");
    fprintf(stderr, "    -o           Output count table.\n");
    fprintf(stderr, "    -q           Map Quality to filter bam.\n");
    fprintf(stderr, "    -no-header   Ignore header in the output.\n");
    return 1;
}

static struct args {
    const char *input_fname;
    const char *output_fname;
    const char *barcode_fname;
    const char *cb_tag; // attribute in BAM
    const char *group_tag;
    const char *list_fname;
    int n_tag;
    char **tags;
    int dedup;
    int qual_thres;
    int ignore_header;
} args = {
    .input_fname   = NULL,
    .output_fname  = NULL,
    .barcode_fname = NULL,
    .cb_tag        = NULL,
    .group_tag     = NULL,
    .n_tag         = 0,
    .tags          = NULL,
    .dedup         = 0,
    .qual_thres    = 0,
    .ignore_header = 0,
};

static int parse_args(int argc, char **argv)
{
    int i;
    const char *tag  = NULL;
    const char *qual = NULL;
    for (i = 1; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;
        if (strcmp(a, "-o") == 0 ) var = &args.output_fname;
        else if (strcmp(a, "-cb") == 0) var = &args.cb_tag;
        else if (strcmp(a, "-tags") == 0) var = &tag;
        else if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return 1;
        else if (strcmp(a, "-dedup") == 0) {
            args.dedup = 1;
            continue;
        }
        else if (strcmp(a, "-group") == 0) var = &args.group_tag;
        else if (strcmp(a, "-list") == 0) var = &args.list_fname;
        else if (strcmp(a, "-q") == 0) var = &qual;
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

    CHECK_EMPTY(args.cb_tag, "-cb must be set.");
    CHECK_EMPTY(args.input_fname, "Input bam must be set.");
    
    if (tag) {
        kstring_t str = {0,0,0};
        kputs(tag, &str);
        int n = 0;
        int *s = ksplit(&str, ',', &n);
        args.tags = malloc(n*sizeof(char*));
        args.n_tag = n;
        int j;
        for (j = 0; j < n; ++j) {
            args.tags[j] = strdup(str.s + s[j]);
            if (strlen(args.tags[j]) != 2) error("Unknown tag format. Only two character allowed.");
        }
        free(s);
        free(str.s);
    }
    
    if (qual) args.qual_thres = str2int((char*)qual);
    if (args.qual_thres) args.qual_thres = 0;
    
    return 0;
}
struct counts {
    int n_group; // group inited of this individual
    struct dict ***counts_per_group; // [n_tag][n_group] point struct dict
};

static void counts_destroy(struct counts *c)
{
    int i;
    for (i = 0; i < args.n_tag; ++i) {
        struct dict **c1 = c->counts_per_group[i];
        int j;
        for (j = 0; j < c->n_group; ++j) {
            if (c1[j] == NULL) continue;
            dict_destroy(c1[j]);
        }
        free(c1);
    }
    free(c->counts_per_group);
}
    
int bam_count_attr(int argc, char *argv[])
{
    double t_real;
    t_real = realtime();

    if (parse_args(argc, argv)) return usage();
    
    htsFile *fp  = hts_open(args.input_fname, "r");
    CHECK_EMPTY(fp, "%s : %s.", args.input_fname, strerror(errno));
    htsFormat type = *hts_get_format(fp);
    if (type.format != bam && type.format != sam)
        error("Unsupported input format, only support BAM/SAM/CRAM format.");
    bam_hdr_t *hdr = sam_hdr_read(fp);
    CHECK_EMPTY(hdr, "Failed to open header.");
    FILE *out = args.output_fname == NULL ? stdout : fopen(args.output_fname, "w");
    CHECK_EMPTY(out, "%s : %s.", args.output_fname, strerror(errno));

    int is_dyn_alloc = 1;

    struct dict *bc_dict = dict_init();
    struct dict *group_dict = NULL;
    if (args.group_tag) group_dict = dict_init();
    
    int n_alloc = 0;
    int m_alloc = 0;
    struct counts *counts = NULL;

    if (args.barcode_fname) {
        dict_read(bc_dict, args.barcode_fname);
        counts = malloc( dict_size(bc_dict) * sizeof(struct counts));
        int i;
        for (i = 0; i < dict_size(bc_dict); ++i) {
            counts[i].n_group = 0;
            counts[i].counts_per_group = calloc(args.n_tag, sizeof(void*));
        }
        is_dyn_alloc = 0;
    }
    
    bam1_t *b;
    bam1_core_t *c;
    int ret;
    b = bam_init1();
    c = &b->core;
    
    while ((ret = sam_read1(fp, hdr, b)) >= 0) {

        if (b->core.qual < args.qual_thres) continue;
        
        uint8_t *tag = bam_aux_get(b, args.cb_tag);
        if (!tag) continue; // skip records without cell Barcodes
        
        char *name = (char*)(tag+1);
        
        int id = -1; // individual index
        
        if (is_dyn_alloc == 0) {
            id = dict_query(bc_dict, name);
            if (id == -1) continue;
        }
        else {
            id = dict_push(bc_dict, name);
            if (id == n_alloc) {
                if (n_alloc == m_alloc) {
                    m_alloc = m_alloc == 0 ? 1024 : m_alloc<<1;
                    counts = realloc(counts, sizeof(struct counts)*m_alloc);
                }
                memset(&counts[id], 0, sizeof(struct counts));
                counts[id].counts_per_group = calloc(args.n_tag, sizeof(void*));
                n_alloc++;
            }
        }
        
        struct counts *cnt = &counts[id];
        
        // dynamic allocate group tag
        int grp_id = 0; // group index
        
        if (args.group_tag) {
            uint8_t *tag = bam_aux_get(b, args.group_tag);
            if (!tag) continue;                
            grp_id = dict_push(group_dict, (char*)(tag+1));
        }
        int new_group = grp_id+1;
        if (cnt->n_group < new_group) {
            int i;
            for (i = 0; i < args.n_tag; ++i) {
                cnt->counts_per_group[i] = realloc(cnt->counts_per_group[i], new_group*sizeof(void *));
                int j;
                for (j = cnt->n_group; j < new_group; ++j)
                    cnt->counts_per_group[i][j] = NULL;
            }
            cnt->n_group = new_group;
        }
        
        int i;
        for (i = 0; i < args.n_tag; ++i) {
            uint8_t *va = bam_aux_get(b, args.tags[i]);
            if (!va) continue;
            struct dict *d = cnt->counts_per_group[i][grp_id];
            if (d == NULL) d = dict_init();
            dict_push(d, (char*)(va+1));
            cnt->counts_per_group[i][grp_id] = d;
        }
    }

    bam_destroy1(b);
    bam_hdr_destroy(hdr);
    sam_close(fp);

    if (ret != -1) warnings("Truncated file?");

    // header
    if (args.ignore_header == 0) {
        fputs("BARCODE\tRaw", out);
        if (args.group_tag) {
            int i;
            for (i = 0; i < dict_size(group_dict); ++i) {
                int j;
                for (j = 0; j < args.n_tag; ++j) {
                    fputc('\t', out);
                    fputs(dict_name(group_dict,i), out);
                    fputc('_', out);
                    fputs(args.tags[j], out);
                }            
            }
        }
        else {
            int i;
            for (i = 0; i < args.n_tag; ++i) {
                fputc('\t', out);
                fputs(args.tags[i], out);
            }
        }
        fputc('\n', out);
    }

    int i;
    for (i = 0; i < dict_size(bc_dict); ++i) {
        struct counts *cnt = &counts[i];
        fprintf(out, "%s\t%u", dict_name(bc_dict, i), dict_count(bc_dict, i) );
        if (args.group_tag) {
            int j;
            for (j = 0; j < dict_size(group_dict); ++j) {
                int k;
                for (k = 0; k < args.n_tag; ++k)
                    fprintf(out, "\t%u",
                            cnt->counts_per_group[j][k] == NULL ? 0 :
                            args.dedup == 1 ? dict_size(cnt->counts_per_group[j][k]) :
                            dict_count_sum(cnt->counts_per_group[j][k]));
            }
        }
        else {
            int j;
            for (j = 0; j < args.n_tag; ++j)
                fprintf(out, "\t%u",
                        cnt->counts_per_group[0][j] == NULL ? 0 :
                        args.dedup == 1 ? dict_size(cnt->counts_per_group[0][j]) :
                        dict_count_sum(cnt->counts_per_group[0][j]));
            
        }
        fputc('\n', out);
        counts_destroy(cnt);
    }
    free(counts);
    
    fclose(out);
    dict_destroy(bc_dict);
    if (args.group_tag) dict_destroy(group_dict);
    LOG_print("Real time: %.3f sec; CPU: %.3f sec", realtime() - t_real, cputime());
    return 0;    
}


        
