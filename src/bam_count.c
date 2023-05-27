// count reads/fragments matrix for single-cell datasets
#include "utils.h"
#include "number.h"
#include "dict.h"
#include "dna_pool.h"
#include "htslib/khash.h"
#include "htslib/kstring.h"
#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include "htslib/thread_pool.h"
#include "htslib/hts_endian.h"
#include "pisa_version.h" // mex output
#include "read_tags.h"
// from v0.10, -ttype supported
#include "read_anno.h"
#include "biostring.h"
#include "bam_pool.h"

// accept file list
#include "bam_files.h"

static struct args {
    const char *input_fname;
    const char *whitelist_fname;
    const char *output_fname;
    const char *outdir; // v0.4, support Market EXchange format (MEX) for sparse matrices
        
    struct dict *tags; // a cell barcode tag or two tags for spatial coordinates
    const char *anno_tag; // feature tag
    const char *umi_tag;

    const char *prefix;

    const char *sample_list;
    
    struct dict *features;
    struct dict *barcodes;
    
    int mapq_thres;
    int use_dup;
    int enable_corr_umi;
    int one_hit;

    int stereoseq;

    int n_thread;
    int chunk_size;

    int bin_size;
    
    //uint64_t n_record;  
    uint64_t n_record1; //  records in spliced matrix
    uint64_t n_record2; //  records in unspliced
    uint64_t n_record3; //  records in spanning
    uint64_t n_record4; //  records in antisense
    //int alias_file_cb;
    
    const char *region_type_tag;
    int n_type;
    enum exon_type *region_types;

    int velocity;

    int antisense;
    
    struct bam_files *files;
} args = {
    .input_fname     = NULL,
    .whitelist_fname = NULL,
    .output_fname    = NULL,
    .outdir          = NULL,
    .tags            = NULL,
    .anno_tag        = NULL,
    .umi_tag         = NULL,

    .prefix          = NULL,
    .barcodes        = NULL,
    .features        = NULL,

    .mapq_thres      = 20,
    .use_dup         = 0,
    .enable_corr_umi = 0,
    .one_hit         = 0,

    .stereoseq       = 0,
    
    .n_thread        = 5,
    .chunk_size      = 1000000,
    .bin_size        = 1,
    
    //.fp_in           = NULL,
    //.hdr             = NULL,
    //.n_record        = 0,
    .n_record1       = 0,
    .n_record2       = 0,
    .n_record3       = 0,
    .n_record4       = 0,
    //.alias_file_cb   = 0,
    
    .region_type_tag = "RE",
    .n_type          = 0,
    .region_types    = NULL,
    .velocity        = 0,
    .antisense       = 0,
    .files           = NULL,
};

struct counts {
    uint32_t count;
    struct PISA_dna_pool *p;

    uint32_t unspliced;
    struct PISA_dna_pool *up;
    
    uint32_t spanning;
    struct PISA_dna_pool *sp;

    uint32_t antisense;
    struct PISA_dna_pool *as;
};

static void memory_release()
{
    //bam_hdr_destroy(args.hdr);
    //sam_close(args.fp_in);
    if (args.tags) dict_destroy(args.tags);
    
    close_bam_files(args.files);
    
    int i;
    int n_feature;
    n_feature = dict_size(args.features);
    for (i = 0; i < n_feature; ++i) {
        struct PISA_dna_pool *v = dict_query_value(args.features, i);
        PISA_idx_destroy(v);
    }
    dict_destroy(args.features);
    dict_destroy(args.barcodes);
}

extern int bam_count_usage();

static int parse_args(int argc, char **argv)
{
    int i;
    const char *mapq = NULL;
    const char *n_thread = NULL;
    const char *region_types = NULL;
    const char *tag_str = NULL;
    const char *bin_size = NULL;
    const char *chunk_size = NULL;

    for (i = 1; i < argc;) {
        const char *a = argv[i++];
        const char **var = 0;
        if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return 1;
        if (strcmp(a, "-tags") == 0 || strcmp(a, "-cb") == 0) var = &tag_str;
        else if (strcmp(a, "-anno-tag") == 0) var = &args.anno_tag;
        else if (strcmp(a, "-list") == 0) var = &args.whitelist_fname;
        else if (strcmp(a, "-umi") == 0) var = &args.umi_tag;
        else if (strcmp(a, "-o") == 0) var = &args.output_fname;
        else if (strcmp(a, "-outdir") == 0) var = &args.outdir;
        else if (strcmp(a, "-q") == 0) var = &mapq;
        else if (strcmp(a, "-t") == 0 || strcmp(a, "-@") == 0) var = &n_thread;
        else if (strcmp(a, "-ttag") == 0) var = &args.region_type_tag;
        else if (strcmp(a, "-ttype") == 0) var = &region_types;
        else if (strcmp(a, "-prefix") == 0) var = &args.prefix;
        else if (strcmp(a, "-sample-list") == 0) var = &args.sample_list;
        else if (strcmp(a, "-bin") == 0) var = &bin_size;
        else if (strcmp(a, "-chunk-size") == 0) var = &chunk_size;
        else if (strcmp(a, "-dup") == 0) {
            args.use_dup = 1;
            continue;
        }
        else if (strcmp(a, "-velo") == 0) {
            args.velocity = 1;
            continue;
        }
        else if (strcmp(a, "-as") == 0) {
            args.antisense = 1;
            continue;
        }
        else if (strcmp(a, "-one-hit") == 0) {
            args.one_hit = 1;
            continue;
        }
        else if (strcmp(a, "-stereoseq") == 0) {
            args.stereoseq = 1;
            continue;
        }
        else if (strcmp(a, "-corr") == 0) {
            //args.enable_corr_umi = 1;
            warnings("Option -corr has been removed since v0.8, to correct UMIs please use `PISA corr` instead.");
            continue;
        }

        else if (strcmp(a, "-file-barcode") ==0) {
            //args.alias_file_cb = 1;
            error("-file-barcode is removed since v0.12.");
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
        error("Unknown argument, %s", a);
    }

    if (args.input_fname == 0 && args.sample_list == NULL) error("No input bam.");
    if (args.input_fname && args.sample_list) error("Input bam conflict with -sample-list.");
    
    if (args.output_fname) {
        warnings("PISA now support MEX format. Old cell X gene expression format is very poor performance. Try -outdir instead of -o.");
    }
    
    if (tag_str == 0) // && args.alias_file_cb == 0)
        error("No cell barcode specified disabled.");

    args.tags = str2tag(tag_str);
    
    if (args.anno_tag == 0) error("No anno tag specified.");

    if (n_thread) args.n_thread = str2int((char*)n_thread);
    if (chunk_size) args.chunk_size = str2int((char*)chunk_size);

    if (args.outdir) {
         struct stat sb;
         if (stat(args.outdir, &sb) != 0) error("Directory %s is not exist.", args.outdir);
         if (S_ISDIR(sb.st_mode) == 0)  error("%s does not look like a directory.", args.outdir);
    }

    if (args.input_fname) {
        args.files = init_bam_line(args.input_fname, args.n_thread > 10 ? 10 : args.n_thread);
    }
    else if (args.sample_list) {
        args.files = init_bam_list(args.sample_list, args.n_thread > 10 ? 10 : args.n_thread);
    }
    else error("Not found input bam file.");

    if (mapq) {
        args.mapq_thres = str2int(mapq);        
    }
    args.features = dict_init();
    dict_set_value(args.features);
    
    args.barcodes = dict_init();

    if (args.whitelist_fname) {
        dict_read(args.barcodes, args.whitelist_fname, 1);
        if (dict_size(args.barcodes) == 0) error("Barcode list is empty?");
    }

    if (region_types) {
        kstring_t str = {0,0,0};
        int n = 0;
        kputs(region_types, &str);
        int *s = str_split(&str, &n);
        if (n == 0) error("Failed to parse -ttype, %s", region_types);
        args.n_type = n;
        args.region_types = malloc(n*sizeof(enum exon_type));
        int k;
        for (k = 0; k < n; ++k) {
            char *rt = str.s+s[k];
            if (strlen(rt) != 1) error("Failed to parse -ttype, %s", region_types);
            enum exon_type type = RE_type_map(rt[0]);
            if (type == type_unknown) error("Unknown type %s", rt);
            args.region_types[k] = type;
        }
        free(s);
        free(str.s);
    }

    if (bin_size) args.bin_size = str2int(bin_size);

    assert(args.bin_size >=1);
    return 0;
}
// decode UMI hex to [ACGT]s
char *stereoseq_decode(char *str, int length)
{
    kstring_t tmp = {0,0,0};
    int l = strlen(str);
    if (l*2 > length) error("Decode string longer than expect.");
    int i;
    for (i = 0; i < l; ++i) {
        switch(str[i]) {
        case '0':
            kputs("AA", &tmp);
            break;
        case '1':
            kputs("CA", &tmp);
            break;
        case '2':
            kputs("GA", &tmp);
            break;
        case '3':
            kputs("TA", &tmp);
            break;
        case '4':
            kputs("AC", &tmp);
            break;
        case '5':
            kputs("CC", &tmp);
            break;
        case '6':
            kputs("GC", &tmp);
            break;
        case '7':
            kputs("TC", &tmp);
            break;
        case '8':
            kputs("AG", &tmp);
            break;
        case '9':
            kputs("CG", &tmp);
            break;
        case 'A':
        case 'a':
            kputs("GG", &tmp);
            break;
        case 'B':
        case 'b':
            kputs("TG", &tmp);
            break;
        case 'C':
        case 'c':
            kputs("AT", &tmp);
            break;
        case 'D':
        case 'd':
            kputs("CT", &tmp);
            break;
        case 'E':
        case 'e':
            kputs("GT", &tmp);
            break;
        case 'F':
        case 'f':
            kputs("TT", &tmp);
            break;
        default:
            error("Invalid hexadecimal digit %c", str[i]);
        }
    }
    
    for ( ; i < length/2; ++i) {
        kputs("AA", &tmp);
    }

    return tmp.s;
}

struct ret {
    struct dict *features;
    struct dict *barcodes;
};
void merge_counts(struct ret *ret)
{
    if (ret == NULL) return;
    int i;
    for (i = 0; i < dict_size(ret->features); ++i) {
        // debug_print("Merging %d features ..", dict_size(ret->features));
        char *feature = dict_name(ret->features, i);
        int idx = dict_query(args.features, feature);
        if (idx < 0) {
            idx = dict_push(args.features, feature);
        }
        struct PISA_dna_pool *v0 = dict_query_value(ret->features, i);
        struct PISA_dna_pool *v = dict_query_value(args.features, idx);
        if (v == NULL) {
            v = PISA_dna_pool_init();
            dict_assign_value(args.features, idx, v);
        }
        
        int j;
        for (j = 0; j < v0->l; ++j) {
            struct PISA_dna *d = &v0->data[j];
            if (d->idx == -1) continue; // cell barcode cached but no record
            //debug_print("%d", d->idx);
            char *barcode = dict_name(ret->barcodes,d->idx);
            assert(barcode);
            int cell_id = dict_query(args.barcodes, barcode);
            if (cell_id == -1) {
                cell_id = dict_push(args.barcodes, barcode);
            }
            struct PISA_dna *c = PISA_idx_query(v, cell_id);
            if (c == NULL) {
                c = PISA_idx_push(v, cell_id);
                struct counts *counts = malloc(sizeof(struct counts));
                memset(counts, 0, sizeof(struct counts));
                counts->p = PISA_dna_pool_init();
                if (args.umi_tag) {
                    counts->up = args.velocity == 1 ? PISA_dna_pool_init() : NULL;
                    counts->sp = args.velocity == 1 ? PISA_dna_pool_init() : NULL;
                    counts->as = args.velocity == 1 ? PISA_dna_pool_init() : NULL;
                }
                c->data = counts;
            }
            
            struct counts *counts = c->data;
            struct counts *c0 = d->data;

            if (args.umi_tag) {
                PISA_pool_merge(counts->p, c0->p);
                
                if (args.velocity == 1) {
                    PISA_pool_merge(counts->up, c0->up);
                    PISA_pool_merge(counts->sp, c0->sp);
                    PISA_pool_merge(counts->as, c0->as);
                    
                    PISA_dna_destroy(c0->up);
                    PISA_dna_destroy(c0->sp);
                    PISA_dna_destroy(c0->as);
                }
                PISA_dna_destroy(c0->p);
            } else {
                counts->count += c0->count;
                if (args.velocity == 1) {
                    counts->unspliced += c0->unspliced;
                    counts->spanning += c0->spanning;
                }
            }
            free(c0);
        }
        PISA_idx_destroy(v0);
    }
    dict_destroy(ret->features);
    dict_destroy(ret->barcodes);
    free(ret);
}
//copy from sam.c
static inline int aux_type2size(uint8_t type)
{
    switch (type) {
    case 'A': case 'c': case 'C':
        return 1;
    case 's': case 'S':
        return 2;
    case 'i': case 'I': case 'f':
        return 4;
    case 'd':
        return 8;
    case 'Z': case 'H': case 'B':
        return type;
    default:
        return 0;
    }
}
static inline uint8_t *skip_aux(uint8_t *s, uint8_t *end)
{
    int size;
    uint32_t n;
    if (s >= end) return end;
    size = aux_type2size(*s); ++s; // skip type
    switch (size) {
    case 'Z':
    case 'H':
        while (s < end && *s) ++s;
        return s < end ? s + 1 : end;
    case 'B':
        if (end - s < 5) return NULL;
        size = aux_type2size(*s); ++s;
        n = le_to_u32(s);
        s += 4;
        if (size == 0 || end - s < size * n) return NULL;
        return s + size * n;
    case 0:
        return NULL;
    default:
        if (end - s < size) return NULL;
        return s + size;
    }
}

static char *retrieve_tags(bam1_t *b, struct dict *tags)
{
    int l = dict_size(tags);

    if (l == 1) {
        uint8_t *data = bam_aux_get(b, dict_name(tags,0));
        if (!data) return NULL;
        return strdup((char*)(data+1));
    }
    
    kstring_t str = {0,0,0};
    char tag[3]; tag[2] = '\0';

    uint8_t *s, *end;
    char **vals = malloc(l*sizeof(char**));
    memset(vals,0, l*sizeof(char**));
    
    s = bam_get_aux(b);
    end = b->data + b->l_data;
    int count = 0;
    while (s != NULL && end - s >= 3) {
        if (count == l) break;
        tag[0] = *s++;
        tag[1] = *s++;
        int idx = dict_query(tags, tag);
        if (idx != -1) {
            // Check the tag value is valid and complete
            uint8_t *e = skip_aux(s, end);
            if ((*s == 'Z' || *s == 'H') && *(e - 1) != '\0') {
                error("Corrupted aux data for read %s", bam_get_qname(b));
            }
            if (e != NULL) {
                // put empty space to unconverted value
                kstring_t tmp = {0,0,0};

                if (*s == 'C' || *s == 'c') {
                    uint8_t va = bam_aux2i(s);
                    if (args.bin_size > 1) va = (uint8_t)(va/args.bin_size)*args.bin_size;
                    kputw(va, &tmp);
                } else if (*s == 'S' || *s == 's') {
                    uint16_t va = bam_aux2i(s);
                    if (args.bin_size > 1) va = (uint16_t)(va/args.bin_size)*args.bin_size;
                    kputw(va, &tmp);
                } else if (*s == 'i' || *s == 'I') {
                    uint32_t va = bam_aux2i(s);
                    if (args.bin_size > 1) va = (uint32_t)(va/args.bin_size)*args.bin_size;
                    kputw(va, &tmp);
                } else if (*s == 'f' || *s == 'd') {
                    double va = bam_aux2f(s);                    
                    kputd(va, &tmp);
                } else if (*s == 'H' || *s == 'Z') {
                    char *va = bam_aux2Z(s);
                    kputs(va, &tmp);
                }
               
                vals[idx] = tmp.s; //strndup((char*)(s+1), s0-s-1);
                count ++;
            } else {
                error("Corrupted aux data for read %s", bam_get_qname(b));
            }
        }
        
        s = skip_aux(s, end);
    }
    
    int i;
    for (i = 0; i < l; ++i) {
        if (vals[i] == NULL) {
            int j;
            for (j = i+1; j < l; ++j) {
                if (vals[j]) free(vals[j]);
            }
            free(vals);
            if (str.m) free(str.s);
            return NULL;
        }
        if (i) kputc('\t', &str);
        kputs(vals[i], &str);
        free(vals[i]);
    }

    free(vals);
    return str.s;
}
//struct ret *count_matrix_core(bam1_t *b, char *tag)
static void *run_it(void *_p)
{
    struct bam_pool *p = (struct bam_pool*)_p;
    if (p == NULL) return NULL;
    struct ret *ret = malloc(sizeof(*ret));
    ret->features = dict_init();
    ret->barcodes = dict_init();
    dict_set_value(ret->features);
    dict_set_value(ret->barcodes);

    if (args.whitelist_fname) {
        dict_read(ret->barcodes, args.whitelist_fname, 1);
    }
    
    kstring_t tmp = {0,0,0};
    kstring_t str = {0,0,0};

    int record;
    for (record = 0; record < p->n; ++record) {
        tmp.l = 0;
        str.l = 0;
        
        bam1_t *b = &p->bam[record];
        if (args.n_type > 0) {
            uint8_t *data = bam_aux_get(b, args.region_type_tag);
            if (!data) continue;
            int region_type_flag = 0;
            int k;
            for (k = 0; k < args.n_type; ++k) {
                if (args.region_types[k] == RE_type_map(data[1])) {
                    region_type_flag = 1;
                    break;
                }
            }
            if (region_type_flag == 0) continue;
        }
                
        uint8_t *anno_tag = bam_aux_get(b, args.anno_tag);
        //if (!anno_tag) goto skip_this_record;
        if (!anno_tag) continue;
        
        if (args.umi_tag) {
            uint8_t *umi_tag = bam_aux_get(b, args.umi_tag);
            // if (!umi_tag) goto skip_this_record;
            if (!umi_tag) continue;
        }

        int unspliced = 0;
        int spanning = 0;
        int antisense = 0; // v0.12
        if (args.velocity) {
            uint8_t *data = bam_aux_get(b, args.region_type_tag);
            // if (!data) goto skip_this_record;
            if (!data) continue;
            if (RE_type_map(data[1]) == type_unknown)  continue; //goto skip_this_record;
            if (RE_type_map(data[1]) == type_ambiguous) continue; //goto skip_this_record;
            if (RE_type_map(data[1]) == type_intergenic) continue; // goto skip_this_record;

            // update 2023/03/20, exonintron labeled for spanning and unspliced matrix
            if (RE_type_map(data[1]) == type_exon_intron) spanning = 1;
            
            if (RE_type_map(data[1]) == type_exon_intron || RE_type_map(data[1]) == type_intron) unspliced = 1;
            else if (RE_type_map(data[1]) == type_exon_intron) spanning = 1;
            else if (RE_type_map(data[1]) == type_antisense) antisense = 1;
            else if (RE_type_map(data[1]) == type_antisense_intron) antisense = 1;
        }
        
        //if (args.tags) {

        // todo: perform improvement
        char *tag = retrieve_tags(b,args.tags);
        /*
        int i;
        for (i = 0; i < dict_size(args.tags); ++i) {
            uint8_t *tag0 = bam_aux_get(b, dict_name(args.tags,i));
            if (!tag0) {
                tmp.l = 0;
                break;
            }
                    
            if (tmp.l) kputc('\t', &tmp);
            if (*tag0 == 'S' || *tag0 == 's' || *tag0 == 'c' || *tag0 == 'i' || *tag0 == 'I') {
                int64_t va = bam_aux2i(tag0);
                kputw(va, &tmp);
            } else if (*tag0 == 'f' || *tag0 == 'd') {
                double va = bam_aux2f(tag0);
                kputd(va, &tmp);
            } else if (*tag0 == 'H' || *tag0 == 'Z') {
                char *va = bam_aux2Z(tag0);
                kputs(va, &tmp);
            }
            
        }
        
        if (tmp.l == 0) continue;
        tag = tmp.s;
        */

        if (tag == NULL) continue;
        
        int cell_id;
        if (args.whitelist_fname) {
            cell_id = dict_query(ret->barcodes, tag);
            free(tag);
            if (cell_id == -1) continue; // goto skip_this_record;
        }
        else {
            cell_id = dict_push(ret->barcodes, tag);
            free(tag);
        }
        
        // for each feature
        kputs((char*)(anno_tag+1), &str);
        int n_gene;
        int *s = str_split(&str, &n_gene); // seperator ; or ,
        
        // Sometime two or more genes or functional regions can overlapped with each other, in default PISA counts the reads for both of these regions.
        // But if -one-hit set, these reads will be filtered.
        if (args.one_hit == 1 && n_gene >1) {
            // free(str.s);
            free(s);
            // goto skip_this_record;
            continue;
        }
        
        int i;
        for (i = 0; i < n_gene; ++i) {
            // Features (Gene or Region)
            char *val = str.s + s[i];
            
            int idx = dict_query(ret->features, val);
            if (idx == -1) idx = dict_push(ret->features, val);

            struct PISA_dna_pool *v = dict_query_value(ret->features, idx);

            if (v == NULL) {
                v = PISA_dna_pool_init();
                dict_assign_value(ret->features, idx, v);
            }
            // not store cell barcode for each hash, use id number instead to reduce memory
            struct PISA_dna *c= PISA_idx_query(v, cell_id);
            if (c == NULL) {
                c = PISA_idx_push(v, cell_id);
                //if (c->data == NULL) {
                struct counts *counts = malloc(sizeof(struct counts));
                memset(counts, 0, sizeof(struct counts));
                
                if (args.umi_tag)  {
                    counts->p = PISA_dna_pool_init();
                    counts->up = args.velocity == 1 ? PISA_dna_pool_init() : NULL;
                    counts->sp = args.velocity == 1 ? PISA_dna_pool_init() : NULL;
                    counts->as = args.velocity == 1 ? PISA_dna_pool_init() : NULL;
                }
                
                c->data = counts;
            }
            
            if (args.umi_tag) {
                uint8_t *umi_tag = bam_aux_get(b, args.umi_tag);
                assert(umi_tag);
                char *val = (char*)(umi_tag+1);
                assert(c->data);
                
                char *val0 = NULL;
                if (args.stereoseq) {
                    val0 = stereoseq_decode(val, 10);
                }
                
                struct counts *count = c->data;
                
                if (args.antisense && antisense) {
                    PISA_dna_push(count->as, val0 ? val0 : val); 
                } else {
                    // total
                    PISA_dna_push(count->p, val0 ? val0 : val);
                    
                    if (args.velocity && unspliced)
                        PISA_dna_push(count->up, val0 ? val0 : val);
                    
                    if (args.velocity && spanning)
                        PISA_dna_push(count->sp, val0 ? val0 : val);
                }
                
                if (val0) free(val0);
            }
            else {
                struct counts *count = c->data;
                count->count++;
                
                if (args.velocity && unspliced)
                    count->unspliced++;
                
                if (args.velocity && spanning)
                    count->spanning++;
                
                if (args.antisense && antisense)
                    count->antisense++;
            }
        }
        free(s);
    }

    if (str.m) free(str.s);
    if (tmp.m) free(tmp.s);

    bam_pool_destory(p);
    
    return ret;
}

static void update_counts()
{
    int n_feature = dict_size(args.features);
    int i;
    for (i = 0; i < n_feature; ++i) {
        struct PISA_dna_pool *v = dict_query_value(args.features, i);
        int j;
        int n_cell = v->l;
        for (j = 0; j < n_cell; ++j) {
            struct counts *count = v->data[j].data;
            assert(count);
            if (args.umi_tag) {
                int size = count->p->l;
                PISA_dna_destroy(count->p);
                count->count = size;
            }

            if (args.velocity) {
                if (args.umi_tag) {
                    int size = count->up->l;
                    PISA_dna_destroy(count->up);
                    count->unspliced = size;

                    size = count->sp->l;
                    PISA_dna_destroy(count->sp);
                    count->spanning = size;

                    size = count->as->l;
                    PISA_dna_destroy(count->as);
                    count->antisense = size;
                }
            }
            // args.n_record += count->count;
            // args.n_record2 += count->unspliced;
            // args.n_record3 += count->spanning;
            // if (count->count > 0 && count->count != count->unspliced) args.n_record1++;
            if (count->count > count->unspliced) args.n_record1++;
            if (count->unspliced > 0) args.n_record2++;
            if (count->spanning > 0) args.n_record3++;
            if (count->antisense > 0) args.n_record4++;
        }
    }
}
static void write_outs()
{
    int n_barcode = dict_size(args.barcodes);
    int n_feature = dict_size(args.features);

    if (n_barcode == 0) error("No barcode found.");
    if (n_feature == 0) error("No feature found.");
    if (args.n_record1 == 0) {
        warnings("No anntated record found.");
        return;
    }

    if (args.outdir) {
        kstring_t barcode_str = {0,0,0};
        kstring_t feature_str = {0,0,0};
        kstring_t mex_str = {0,0,0};
        kstring_t unspliced_str = {0,0,0};
        kstring_t spanning_str = {0,0,0};
        kstring_t antisense_str = {0,0,0};
    
        kputs(args.outdir, &barcode_str);
        kputs(args.outdir, &feature_str);
        kputs(args.outdir, &mex_str);
        kputs(args.outdir, &unspliced_str);
        kputs(args.outdir, &spanning_str);
        kputs(args.outdir, &antisense_str);
        
        if (args.outdir[strlen(args.outdir)-1] != '/') {
            kputc('/', &barcode_str);
            kputc('/', &feature_str);
            kputc('/', &mex_str);
            kputc('/', &unspliced_str);
            kputc('/', &spanning_str);
            kputc('/', &antisense_str);
        }
    
        if (args.prefix) {
            kputs(args.prefix, &barcode_str);
            kputs(args.prefix, &feature_str);
            kputs(args.prefix, &mex_str);
            kputs(args.prefix, &unspliced_str);
            kputs(args.prefix, &spanning_str);
            kputs(args.prefix, &antisense_str);
        }
        
        kputs("barcodes.tsv.gz", &barcode_str);
        kputs("features.tsv.gz", &feature_str);
        if (args.velocity)
            kputs("spliced.mtx.gz", &mex_str);
        else
            kputs("matrix.mtx.gz", &mex_str);

        kputs("unspliced.mtx.gz", &unspliced_str);
        kputs("spanning.mtx.gz", &spanning_str);
        kputs("antisense.mtx.gz", &antisense_str);
        
        BGZF *barcode_fp = bgzf_open(barcode_str.s, "w");
        bgzf_mt(barcode_fp, args.n_thread, 256);
        CHECK_EMPTY(barcode_fp, "%s : %s.", barcode_str.s, strerror(errno));
        
        int i;

        kstring_t str = {0,0,0};
        kstring_t str2 = {0,0,0};
        kstring_t str3 = {0,0,0};
        kstring_t str4 = {0,0,0};
        int l;
        if (1) {
            for (i = 0; i < n_barcode; ++i) {
                char *name = dict_name(args.barcodes, i);
                kputs(name, &str);
                /*
                debug_print("%s", name);

                char *s0 = name;
                char *end = name + strlen(name);

                for (; *s0 && s0 != end;) {
                    if (*s0 == ' ') {
                        s0++;
                        uint8_t type = *s0++;
                        if (type == 'c') {
                            int8_t va = le_to_i8(s0);
                            s0 = s0 + 1;
                        } else if (type == 'C') {
                            uint8_t va = *s0;
                            s0 = s0 + 1;
                        } else if (type == 's') {
                            int16_t va = le_to_i16(s0);
                            s0 = s0 + 2;
                        } else if (type == 'S') {
                            uint16_t va = le_to_u16(s0);
                            s0 = s0 + 2;
                        } else if (type == 'i') {
                            int32_t va = le_to_i32(s0);
                            kputw(va, &str);
                            s0 = s0 + 4;
                        } else if (type == 'I') {
                            uint32_t va = le_to_u32(s0);
                            kputw(va, &str);
                            s0 = s0 + 4;
                        } else if (type == 'f') {
                            float va = le_to_float(s0);
                            kputd(va, &str);
                            s0 = s0 + 4;
                        } else if (type == 'd') {
                            double va = le_to_double(s0);
                            kputd(va, &str);
                            s0 = s0 + 8;
                        } else {
                            error("Unknown type : %s", s0);
                        }
                    } else if (*s0 == '\t') {
                        kputc('\t', &str);
                        s0++;
                    } else {
                        char *val = s0;
                        for (;*s0 && s0 != end && *s0 != '\t'; ++s0);
                        kputsn(val, s0- val,&str);
                        kputs("", &str);
                    }
                    }
                */
                kputc('\n', &str);
            }
            l = bgzf_write(barcode_fp, str.s, str.l);
            if (l != str.l) error("Failed to write.");
            bgzf_close(barcode_fp);
        }
        str.l = 0;
        BGZF *feature_fp = bgzf_open(feature_str.s, "w");
        bgzf_mt(feature_fp, args.n_thread, 256);
        CHECK_EMPTY(feature_fp, "%s : %s.", feature_str.s, strerror(errno));
        for (i = 0; i < n_feature; ++i) {
            kputs(dict_name(args.features,i), &str);
            kputc('\n', &str);
        }
        l = bgzf_write(feature_fp, str.s, str.l);
        if (l != str.l) error("Failed to write.");
        
        bgzf_close(feature_fp);

        str.l = 0;

        BGZF *mex_fp = bgzf_open(mex_str.s, "w");
        CHECK_EMPTY(mex_fp, "%s : %s.", mex_str.s, strerror(errno));
        
        bgzf_mt(mex_fp, args.n_thread, 256);
        kputs("%%MatrixMarket matrix coordinate integer general\n", &str);
        kputs("% Generated by PISA ", &str);
        kputs(PISA_VERSION, &str);
        kputc('\n', &str);
        ksprintf(&str, "%d\t%d\t%" PRIu64 "\n", n_feature, n_barcode, args.n_record1);

        BGZF *unspliced_fp = NULL;
        BGZF *spanning_fp = NULL;
        BGZF *antisense_fp = NULL;
        if (args.velocity) {
            unspliced_fp = bgzf_open(unspliced_str.s, "w");
            if (unspliced_fp == NULL) error("%s : %s.", unspliced_str.s, strerror(errno));
        
            bgzf_mt(unspliced_fp, args.n_thread, 256);
            kputs("%%MatrixMarket matrix coordinate integer general\n", &str2);
            kputs("% Generated by PISA ", &str2);
            kputs(PISA_VERSION, &str2);
            kputc('\n', &str2);
            ksprintf(&str2, "%d\t%d\t%" PRIu64 "\n", n_feature, n_barcode, args.n_record2);

            spanning_fp = bgzf_open(spanning_str.s, "w");
            if (spanning_fp == NULL) error("%s : %s.", spanning_str.s, strerror(errno));
        
            bgzf_mt(spanning_fp, args.n_thread, 256);
            kputs("%%MatrixMarket matrix coordinate integer general\n", &str3);
            kputs("% Generated by PISA ", &str3);
            kputs(PISA_VERSION, &str3);
            kputc('\n', &str3);
            ksprintf(&str3, "%d\t%d\t%" PRIu64 "\n", n_feature, n_barcode, args.n_record3);
        }
        if (args.antisense) {
            antisense_fp = bgzf_open(antisense_str.s, "w");
            if (antisense_fp == NULL) error("%s : %s.", antisense_str.s, strerror(errno));
            
            bgzf_mt(antisense_fp, args.n_thread, 256);
            kputs("%%MatrixMarket matrix coordinate integer general\n", &str4);
            kputs("% Generated by PISA ", &str4);
            kputs(PISA_VERSION, &str4);
            kputc('\n', &str4);
            ksprintf(&str4, "%d\t%d\t%" PRIu64 "\n", n_feature, n_barcode, args.n_record4);
        }
        
        for (i = 0; i < n_feature; ++i) {
            struct PISA_dna_pool *v = dict_query_value(args.features, i);
            int j;
            int n_cell = v->l;
            for (j = 0; j < n_cell; ++j) {
                struct counts *count = v->data[j].data;
                if (args.velocity) {
                    int spliced = count->count - count->unspliced;
                    if (spliced > 0) ksprintf(&str, "%d\t%d\t%u\n", i+1, v->data[j].idx+1, spliced);
                    if (count->unspliced > 0) ksprintf(&str2, "%d\t%d\t%u\n", i+1, v->data[j].idx+1, count->unspliced);
                    if (count->spanning > 0)  ksprintf(&str3, "%d\t%d\t%u\n", i+1, v->data[j].idx+1, count->spanning);
                    if (args.antisense && count->antisense > 0) ksprintf(&str4, "%d\t%d\t%u\n", i+1, v->data[j].idx+1, count->antisense);
                }
                else
                    ksprintf(&str, "%d\t%d\t%u\n", i+1, v->data[j].idx+1, count->count);
                
                free(count);
            }

            if (str.l > 100000000) {
                int l = bgzf_write(mex_fp, str.s, str.l);
                if (l != str.l) error("Failed to write file.");
                str.l = 0;

                if (args.velocity) {
                    l = bgzf_write(unspliced_fp, str2.s, str2.l);
                    if (l != str2.l) error("Failed to write file.");
                    str2.l = 0;

                    l = bgzf_write(spanning_fp, str3.s, str3.l);
                    if (l != str3.l) error("Failed to write file.");
                    str3.l = 0;
                }
                if (args.antisense) {
                    l = bgzf_write(antisense_fp, str4.s, str4.l);
                    if (l != str4.l) error("Failed to write file.");
                    str4.l = 0;
                }
            }
        }

        if (str.l) {
            l = bgzf_write(mex_fp, str.s, str.l);
            if (l != str.l) error("Failed to wirte.");
        }

        if (str2.l) {
            l = bgzf_write(unspliced_fp, str2.s, str2.l);
            if (l != str2.l) error("Failed to wirte.");
        }

        if (str3.l) {
            l = bgzf_write(spanning_fp, str3.s, str3.l);
            if (l != str3.l) error("Failed to wirte.");
        }

        if (str4.l) {
            l = bgzf_write(antisense_fp, str4.s, str4.l);
            if (l != str4.l) error("Failed to wirte.");
        }

        free(str.s);
        if (str2.m) free(str2.s);
        if (str3.m) free(str3.s);
        if (str4.m) free(str4.s);
        free(mex_str.s);
        free(barcode_str.s);
        free(feature_str.s);
        if (unspliced_str.m) free(unspliced_str.s);
        if (spanning_str.m) free(spanning_str.s);
        if (antisense_str.m) free(antisense_str.s);
        bgzf_close(mex_fp);
        if (unspliced_fp) bgzf_close(unspliced_fp);
        if (spanning_fp) bgzf_close(spanning_fp);
        if (antisense_fp) bgzf_close(antisense_fp);
    }
    
    // header
    if (args.output_fname) {
        int i;
        FILE *out = fopen(args.output_fname, "w");
        CHECK_EMPTY(out, "%s : %s.", args.output_fname, strerror(errno));
        fputs("ID", out);
        
        for (i = 0; i < n_barcode; ++i)
            fprintf(out, "\t%s", dict_name(args.barcodes, i));
        fprintf(out, "\n");
        uint32_t *temp = malloc(n_barcode*sizeof(int));        
        for (i = 0; i < n_feature; ++i) {
            struct PISA_dna_pool *v = dict_query_value(args.features, i);
            int j;
            int n_cell = v->l;
            memset(temp, 0, sizeof(int)*n_barcode);
            fputs(dict_name(args.features, i), out);
            for (j = 0; j < n_cell; ++j) {
                int idx = v->data[j].idx;
                temp[idx] = v->data[j].count;
            }

            for (j = 0; j < n_barcode; ++j)
                fprintf(out, "\t%u", temp[j]);
            fputc('\n', out);
        }
        fclose(out);
        free(temp);
    }
}

struct bam_pool *read_files_pool(struct bam_files *files, int size)
{
    struct bam_pool *p = bam_pool_init(size);

    int ret;
    
    for (p->n = 0; p->n < p->m;) {
        ret = read_bam_files(files, &p->bam[p->n]);
        if (ret < 0) break;
        bam_hdr_t *hdr = get_hdr(args.files);
        bam1_t *b =  &p->bam[p->n];
        // char *alias = get_alias(args.files);
        // if (args.alias_file_cb == 1 && !alias)
        // error("No alias found for %s", get_fname(args.files));
        
        bam1_core_t *c;
        c = &b->core;

        if (c->tid <= -1 || c->tid > hdr->n_targets || (c->flag & BAM_FUNMAP)) continue;
        if (c->qual < args.mapq_thres) continue;
        if (args.use_dup == 0 && c->flag & BAM_FDUP) continue;
        
        p->n++;
    }

    if (p->n == 0) {
        bam_pool_destory(p);
        return NULL;
    }
    return p;
}

int count_matrix(int argc, char **argv)
{
    double t_real;
    t_real = realtime();
    if (parse_args(argc, argv)) return bam_count_usage();

    hts_tpool *p = hts_tpool_init(args.n_thread);
    hts_tpool_process *q = hts_tpool_process_init(p, args.n_thread*2, 0);
    hts_tpool_result *r;
    
    for (;;) {
        struct bam_pool *pool = read_files_pool(args.files, args.chunk_size);
        if (pool == NULL) break;
        int block;
        do {
            block = hts_tpool_dispatch2(p, q, run_it, pool, 1);
            if ((r = hts_tpool_next_result(q))) {
                struct ret *ret = (struct ret*)hts_tpool_result_data(r);
                //write_out(d);
                merge_counts(ret);
                hts_tpool_delete_result(r, 0);
            }
        }
        while (block == -1);
    }
    
    hts_tpool_process_flush(q);
 
    while ((r = hts_tpool_next_result(q))) {
        struct ret *ret = (struct ret*)hts_tpool_result_data(r);
        //write_out(d);
        merge_counts(ret);
        hts_tpool_delete_result(r, 0);
    }
    hts_tpool_process_destroy(q);
    hts_tpool_destroy(p);

    update_counts();

    write_outs();
    
    memory_release();
    
    LOG_print("Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB.", realtime() - t_real, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
    return 0;
}
int count_matrix1(int argc, char **argv)
{
    double t_real;
    t_real = realtime();
    if (parse_args(argc, argv)) return bam_count_usage();

#pragma omp parallel num_threads(args.n_thread)
    for (;;) {
        
        struct bam_pool *pool = NULL;

#pragma omp critical (read)
        pool = read_files_pool(args.files, args.chunk_size);
        if (pool == NULL) break;

        struct ret *ret = run_it(pool);
        
#pragma omp critical (merge)
        merge_counts(ret);

    }
    
    update_counts();

    write_outs();
    
    memory_release();
    
    LOG_print("Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB.", realtime() - t_real, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
    return 0;
}
