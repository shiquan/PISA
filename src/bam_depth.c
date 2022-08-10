#include "utils.h"
#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "dict.h"
#include "number.h"
#include "bed.h"

static struct args {
    const char *input_fname;
    const char *output_fname;
    const char *bc_list;
    const char *tag; // cell barcode tag
    const char *umi_tag; // umi tag
    const char *region_fname;
    
    struct dict *barcodes;
    htsFile *fp;
    hts_idx_t *idx;
    FILE *out;
    
    struct bed_spec *B;
    int tid;
    int start;
    int end;
    int strand;

    bam_hdr_t *hdr;
    
    int mapq_thres;
    int n_thread;
} args = {
    .input_fname = NULL,
    .output_fname = NULL,
    .bc_list     = NULL,
    .tag         = "CB",
    .umi_tag     = NULL,
    .barcodes    = NULL,
    .region_fname = NULL,
    .fp          = NULL,
    .idx         = NULL,
    .out         = NULL,
    .B           = NULL,
    .tid         = -1,
    .start       = 0,
    .end         = 0,
    .strand      = BED_STRAND_UNK,
    .mapq_thres  = 20,
    .n_thread    = 4,
};

static int parse_args(int argc, char **argv)
{
    int i;
    const char *mapq = NULL;
    const char *region = NULL;
    const char *threads = NULL;
    for (i = 1; i < argc;) {
        const char *a = argv[i++];
        const char **var = 0;
        if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return 1;
        if (strcmp(a, "-tag") == 0) var = &args.tag;
        else if (strcmp(a, "-list") == 0) var = &args.bc_list;
        else if (strcmp(a, "-umi") == 0) var = &args.umi_tag;
        else if (strcmp(a, "-o") == 0) var = &args.output_fname;
        else if (strcmp(a, "-q") == 0) var = &mapq;
        else if (strcmp(a, "-@") == 0) var = &threads;
        else if (strcmp(a, "-bed") == 0) var = &args.region_fname;
        
        if (var != 0) {
            if (i == argc) error("Miss an argument after %s.", a);
            *var = argv[i++];
            continue;
        }

        if (args.input_fname == NULL) {
            args.input_fname = a;
            continue;
        }

        if (region == NULL) {
            region = a;
            continue;
        }
        
        error("Unknown argument, %s", a);
    }

    if (args.input_fname == 0) error("No input bam.");

    if (args.region_fname && region) error("-bed is conflict with region.");
    if (args.region_fname == NULL && region == NULL) error("Required a bed file or region.");
    if (threads) args.n_thread = str2int((char*)threads);
    
    if (mapq) {
        args.mapq_thres = str2int(mapq); 
    }

    if (args.bc_list) {
        args.barcodes = dict_init();
        dict_read(args.barcodes, args.bc_list, 0);
        if (dict_size(args.barcodes) == 0) error("Barcode list is empty?");
    }
    
    args.fp = hts_open(args.input_fname, "r");
    if (args.fp == NULL) error("%s : %s.", args.input_fname, strerror(errno));
    
    hts_set_threads(args.fp, args.n_thread);

    args.idx = sam_index_load(args.fp, args.input_fname);
    if (args.idx == NULL) error("Failed to load bam index of %s", args.input_fname);

    args.hdr = sam_hdr_read(args.fp);
    if (args.output_fname) {
        args.out = fopen(args.output_fname, "w");
        if (args.out == NULL) error("%s : %s.", args.output_fname, strerror(errno));        
    }
    else args.out = stdout;

    // regions
    if (args.region_fname) {
        args.B = bed_read(args.region_fname);
    }
    else {
        kstring_t str = {0,0,0};
        kputs(region, &str);
        int j;
        char *s = str.s;
        for (j = 0; j < str.l; ++j) {
            if (str.s[j] == ':') {
                str.s[j]=0;
                if (args.tid == -1) {
                    args.tid = sam_hdr_name2tid(args.hdr, s);
                    if (args.tid == -1) error("Not found this chromosome at header. %s", str.s);
                    s = str.s+j+1;
                }
                else {
                    args.end = str2int(s);
                    s = str.s + j+1;
                    if (*s == '+') args.strand = BED_STRAND_FWD;
                    else if (*s == '-') args.strand = BED_STRAND_REV;
                    else error("Unknow strand. %s", s) ;
                    break;
                }

            }
            else if (str.s[j] == '-') {
                str.s[j] = 0;
                args.start = str2int(s) -1;
                s = str.s+j+1;
                if (s == NULL) {
                    args.end = args.start+1;
                    args.strand = BED_STRAND_UNK;
                    break;
                }
            }         
        }

        if (args.end == 0 && s) args.end = str2int(s);
        free(str.s);
    }
    
    return 0;
}

static void memory_release()
{
    if (args.barcodes) dict_destroy(args.barcodes);
    
    bam_hdr_destroy(args.hdr);
    hts_idx_destroy(args.idx);
    hts_close(args.fp);
    bed_spec_destroy(args.B);
    
    if (args.out != stdout) fclose(args.out);    
}

struct depth {
    int pos;
    int depth;
    int strand;
    struct dict *bc;
    struct depth *next;
    struct depth *before;
};

static struct depth *depth_init()
{
    struct depth *d = malloc(sizeof(*d));
    memset(d, 0, sizeof(*d));
    return d;
}

int bam2depth(const hts_idx_t *idx, const int tid, const int start, const int end, const int strand)
{
    hts_itr_t *itr = sam_itr_queryi(idx, tid, start, end);
    bam1_t *b = bam_init1();

    int result;

    struct depth *head = NULL;
    struct depth *tail = NULL;

    int strand0;
    
    if (itr) {
        while ((result = sam_itr_multi_next(args.fp, itr, b)) >= 0) {
            bam1_core_t *c = &b->core;
            if (c->qual < args.mapq_thres) continue;
            if (c->flag & BAM_FSECONDARY) continue;
            if (c->flag & BAM_FUNMAP) continue;
            if (c->flag & BAM_FQCFAIL) continue;
            if (c->flag & BAM_FDUP) continue;

            if (c->flag & BAM_FREVERSE) strand0 = BED_STRAND_REV;
            else strand0 = BED_STRAND_FWD;
            
            if (strand == BED_STRAND_FWD && strand0 == BED_STRAND_REV) continue;
            if (strand == BED_STRAND_REV && strand0 == BED_STRAND_FWD) continue;
            
            uint8_t *data = NULL;
            if (args.barcodes) {
                data = bam_aux_get(b, args.tag);
                if (data == NULL) continue;
                int idx = dict_query(args.barcodes, (char*)(data+1));
                if (idx == -1) continue;
            }

            if (args.umi_tag) {
                data = bam_aux_get(b, args.umi_tag);
                if (data == NULL) continue;                
            }

            int i;
            int pos = c->pos+1;
            while (head && head->pos < pos) {
                if (head->bc && head->depth == 0) {
                    head->depth = dict_size(head->bc);
                    dict_destroy(head->bc);
                }
                if (head->pos > start && head->pos <= end)
                    fprintf(args.out, "%s\t%d\t%c\t%d\n",
                            args.hdr->target_name[tid],
                            head->pos,
                            ".+-"[head->strand+1],
                            head->depth);
                struct depth *tmp = head;
                head = head->next;
                
                if (head) head->before = NULL;
                free(tmp);
            }

            for (i = 0; i < b->core.n_cigar; ++i) {
                int cig = bam_cigar_op(bam_get_cigar(b)[i]);
                int ncig = bam_cigar_oplen(bam_get_cigar(b)[i]);
                if (cig == BAM_CMATCH || cig == BAM_CEQUAL || cig == BAM_CDIFF) {
                    int j;
                    for (j = 0; j < ncig; ++j) {
                        if (pos <= start) { pos++; continue;}
                        if (pos > end) { pos++; continue;} // todo:
                        if (head == NULL) {
                            head = depth_init();
                            head->pos = pos;
                            head->next = NULL;
                            head->before = NULL;
                            head->strand = strand0;
                            tail = head;
                        }

                        struct depth *cur;
                        
                        for (cur = tail; cur != NULL; cur = cur->before) {
                            if (cur->pos == pos && cur->strand == strand0) {
                                if (args.umi_tag) {
                                    if (cur->bc == NULL) cur->bc = dict_init();
                                    dict_push(cur->bc, (char*)(data+1));
                                }
                                else cur->depth++;
                                
                                break;
                            }
                            else if (cur->pos < pos) {
                                struct depth *new = depth_init();
                                new->pos = pos;
                                new->next = cur->next;
                                new->before = cur;
                                new->strand = strand0;
                                if (cur->next) cur->next->before = new;
                                cur->next = new;

                                if (args.umi_tag) {
                                    if (new->bc == NULL) new->bc = dict_init();
                                    dict_push(new->bc, (char*)(data+1));
                                }
                                else new->depth++;

                                if (cur == tail) tail = new;
                                
                                break;
                            }   
                        }

                        if (cur == NULL) {
                            struct depth *new = depth_init();
                            new->pos = pos;
                            new->next = head;
                            new->before = NULL;
                            new->strand = strand0;

                            head->before = new;

                            if (args.umi_tag) {
                                if (new->bc == NULL) new->bc = dict_init();
                                dict_push(new->bc, (char*)(data+1));
                            }
                            else new->depth++;
                            head = new;
                        }
                        pos++;
                    }
                }
                else if (cig == BAM_CDEL) {
                    pos += ncig;
                }
                else if (cig == BAM_CREF_SKIP) {
                    pos += ncig;
                }
            }

        }
        
        struct depth *cur = head;
        while (cur) {
            //for (cur = head; cur != NULL; cur = cur->next) {
            if (cur->bc && cur->depth == 0) {
                cur->depth = dict_size(cur->bc);
                dict_destroy(cur->bc);
            }
            if (cur->pos > start && cur->pos <= end)
                fprintf(args.out, "%s\t%d\t%c\t%d\n", args.hdr->target_name[tid], cur->pos, ".+-"[cur->strand+1], cur->depth);
            struct depth *tmp = cur;
            cur = cur->next;
            if (cur) cur->before = NULL;
            free(tmp);
        }
        
        if (result < -1)
            warnings("Failed to retrieve region due to truncated file or corrupt bam index.");
        
        hts_itr_destroy(itr);
    }
    bam_destroy1(b);
    return 0;
}

extern int depth_usage();

int depth_main(int argc, char **argv)
{
    if (parse_args(argc, argv) == 1) return depth_usage();

    if (args.B == NULL && args.tid != -1)
        return bam2depth(args.idx, args.tid, args.start, args.end, args.strand);

    int i;
    for (i = 0; i < args.B->n; ++i) {
        struct bed *bed = &args.B->bed[i];
        char *chr = dict_name(args.B->seqname, bed->seqname);
        int tid = sam_hdr_name2tid(args.hdr, chr);
        bam2depth(args.idx, tid, bed->start, bed->end, bed->strand);
    }

    memory_release();
    
    return 0;
}
