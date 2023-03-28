#include "utils.h"
#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "dict.h"
#include "number.h"
#include "bed.h"
#include "coverage.h"

const static int ws1m = 1000000;

static struct args {
    const char *input_fname;
    const char *output_fname;
    const char *bc_list;
    const char *tag; // cell barcode tag
    const char *umi_tag; // umi tag
    const char *region_fname;
    
    struct dict *barcodes;
    int fix_barcodes;
    
    htsFile *fp;
    hts_idx_t *idx;
    FILE *out; // 

    struct bed_spec *B;
    int tid;
    int start;
    int end;
    int strand;

    bam_hdr_t *hdr;
    
    int mapq_thres;
    int n_thread;

    int print_0;
    int ignore_strand;
    int split_by_tag;
    int alias_tag;
    struct dict *alias;
    int *alias_idx;

    // int bg; // bedGraph like format, chrom, start (0 based), end, depth, strand
} args = {
    .input_fname = NULL,
    .output_fname = NULL,
    .bc_list     = NULL,
    .tag         = NULL,
    .umi_tag     = NULL,
    .barcodes    = NULL,
    .fix_barcodes= 0,
    
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
    .print_0     = 0,
    .ignore_strand = 0,
    .split_by_tag= 0,
    .alias_tag   = 0,
    .alias       = NULL,
    .alias_idx   = NULL,
    // .bg          = 0
};

static void print_node(struct depth *d, int tid, FILE *out)
{
    if (args.ignore_strand) {
        fprintf(out, "%s\t%d\t.\t%d\n", args.hdr->target_name[tid], d->pos, d->dep1);
    } else {
        fprintf(out, "%s\t%d\t+\t%d\n", args.hdr->target_name[tid], d->pos, d->dep1);
        fprintf(out, "%s\t%d\t-\t%d\n", args.hdr->target_name[tid], d->pos, d->dep2);
    }
}
static void print_node_id(struct depth *d, int tid, int id, FILE *out)
{
    char *v;
    if (args.alias_tag) {
        if (id == -1) v = "others";
        else v=dict_name(args.alias, id);
    } else {
        v= dict_name(args.barcodes, id);
    }
    
    if (args.ignore_strand) {
        fprintf(out, "%s\t%d\t.\t%d\t%s\n",args.hdr->target_name[tid], d->pos, d->dep1, v);
    } else {
        fprintf(out, "%s\t%d\t+\t%d\t%s\n",args.hdr->target_name[tid], d->pos, d->dep1, v);
        fprintf(out, "%s\t%d\t-\t%d\t%s\n",args.hdr->target_name[tid], d->pos, d->dep2, v);
    }
}

static void print_node_id0(int pos, int tid, int id, FILE *out)
{
    char *v;
    if (args.alias_tag) {
        if (id == -1) v = "others";
        else v = dict_name(args.alias, id);
    } else {
        v = dict_name(args.barcodes, id);
    }

    if (args.ignore_strand) {
        fprintf(out, "%s\t%d\t.\t0\t%s\n",args.hdr->target_name[tid], pos, v);
    } else {
        fprintf(out, "%s\t%d\t+\t0\t%s\n",args.hdr->target_name[tid], pos, v);
        fprintf(out, "%s\t%d\t-\t0\t%s\n",args.hdr->target_name[tid], pos, v);
    }
}
static void print_node0(int pos, int tid, FILE *out)
{
    if (args.split_by_tag == 0) {
        if (args.ignore_strand) {
            fprintf(out, "%s\t%d\t.\t0\n", args.hdr->target_name[tid], pos);
        } else {
            fprintf(out, "%s\t%d\t+\t0\n", args.hdr->target_name[tid], pos);
            fprintf(out, "%s\t%d\t-\t0\n", args.hdr->target_name[tid], pos);
        }
    } else {
        int id;
        int n = 0;
        if (args.alias_tag) n = dict_size(args.alias);
        else n = dict_size(args.barcodes);
        for (id = 0; id < n; ++id) {
            print_node_id0(pos, tid, id, out);
        }
    }
}

int depth2file(struct depth *d, int tid, int start, int end, FILE *out)
{
    int last = start + 1;
    // struct depth *d = _d;
    while (d) {
        int last_pos = d->pos;
        int last_id = 0;
        while (d && d->pos == last_pos) {
            if (d->pos > start && d->pos <= end) {
                if (args.print_0) {
                    for (; last < d->pos; ++last) {
                        print_node0(last, tid, out);
                    }
                    last = d->pos +1;
                }
                
                if (args.split_by_tag == 0) {
                    print_node(d, tid, out);
                } else {
                    if (args.print_0 && last_id < d->id) {
                        for (; last_id < d->id; last_id++) {
                            print_node_id0(last_pos, tid, last_id, out);
                        }
                    }
                    
                    print_node_id(d, tid, d->id, out);
                    
                last_id = d->id+1;
                }
            }
            
            struct depth *tmp = d;
            d = d->next;
            //d->before = NULL;
            free(tmp);
        }

        if (args.split_by_tag) {
            int n = 0;
            if (args.alias_tag) n = dict_size(args.alias);
            else n = dict_size(args.barcodes);
            for (; args.print_0 && last_id < n; last_id++) {
                print_node_id0(last_pos, tid, last_id, out);
            }
        }
    }

    for (; args.print_0 && last <= end; ++last) {
        print_node0(last, tid, out);
    }
    
    return 0;
}

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
        else if (strcmp(a, "-0") == 0) {
            args.print_0 = 1;
            continue;
        }
        else if (strcmp(a, "-split") == 0) {
            args.split_by_tag = 1;
            continue;
        }
        else if (strcmp(a, "-is") == 0) {
            args.ignore_strand = 1;
            continue;
        }
        /* else if (strcmp(a, "-bg") == 0) { */
        /*     args.bg = 1; */
        /*     continue; */
        /* } */
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
    // if (args.region_fname == NULL && region == NULL) error("Required a bed file or region.");
    
    if (threads) args.n_thread = str2int((char*)threads);
    
    if (mapq) {
        args.mapq_thres = str2int(mapq); 
    }

    if (args.split_by_tag == 1 && args.tag == NULL) error("-tag is required which set -split.");
    if (args.tag) args.barcodes = dict_init();
    
    if (args.bc_list) {
        // args.barcodes = dict_init();
        int val = 0;
        dict_read2(args.barcodes, args.bc_list, &val);
        if (dict_size(args.barcodes) == 0) error("Barcode list is empty?");

        if (args.split_by_tag && val == 1) {
            args.alias_tag = 1;
            int n = dict_size(args.barcodes);
            args.alias = dict_init();
            args.alias_idx = malloc(sizeof(int)*n);

            int i;
            for (i = 0; i < n; ++i) {
                args.alias_idx[i] = -1;
                char *v = dict_query_value(args.barcodes, i);                
                if (v) {
                    int idx = dict_push(args.alias, v);
                    args.alias_idx[i] = idx;
                    assert(idx < dict_size(args.alias));
                    free(v);
                }
            }
        }
        args.fix_barcodes = 1;
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
    if (args.region_fname) { // -bed
        args.B = bed_read(args.region_fname);
        int i;
        for (i = 0; i < args.B->n; ++i) {
            struct bed *bed = &args.B->bed[i];
            char *chr = dict_name(args.B->seqname, bed->seqname);
            int tid = sam_hdr_name2tid(args.hdr, chr);
            struct depth *d = bam2depth(args.idx, tid, bed->start, bed->end, bed->strand,
                                        args.fp, args.mapq_thres, args.ignore_strand,
                                        args.barcodes, args.tag, args.umi_tag, args.split_by_tag,
                                        args.alias_tag, args.alias_idx, args.fix_barcodes);
            depth2file(d, tid, bed->start, bed->end, args.out);
        }
    }
    else if (region) { // target region
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

        if (args.tid == -1 && args.start == 0) {
            args.tid = sam_hdr_name2tid(args.hdr, s);
            if (args.tid == -1) error("Not found this chromosome at header. %s", str.s);
            args.start = 0;
            args.end = args.hdr->target_len[args.tid];
        }
        else if (args.end == 0 && s) args.end = str2int(s);
        
        free(str.s);
        struct depth *d = bam2depth(args.idx, args.tid, args.start, args.end, args.strand,
                                    args.fp, args.mapq_thres, args.ignore_strand,
                                    args.barcodes, args.tag, args.umi_tag, args.split_by_tag,
                                    args.alias_tag, args.alias_idx, args.fix_barcodes);
        depth2file(d, args.tid, args.start, args.end, args.out);
        // depth_destroy(d);        
    }
    else { // whole genome
        int i;
        for (i = 0; i < args.hdr->n_targets; ++i) {
            int len = args.hdr->target_len[i];
            int last = 0;
            int end = last + ws1m;
            if (end > len) end = len;
            
            for (;;) {
                struct depth *d = bam2depth(args.idx, i, last, end, args.strand,
                                            args.fp, args.mapq_thres, args.ignore_strand,
                                            args.barcodes, args.tag, args.umi_tag,
                                            args.split_by_tag, args.alias_tag, args.alias_idx,
                                            args.fix_barcodes);

                depth2file(d, i, last, end, args.out);
                // depth_destroy(d);
                if (end == len) break;
                
                last = end;
                end = end + ws1m;
                if (end > len) end = len;
            }
        }
    }
    
    return 0;
}

static void memory_release()
{
    if (args.barcodes) dict_destroy(args.barcodes);
    
    bam_hdr_destroy(args.hdr);
    hts_idx_destroy(args.idx);
    hts_close(args.fp);
    if (args.B) bed_spec_destroy(args.B);
    
    if (args.out != stdout) fclose(args.out);    
}

extern int depth_usage();

int depth_main(int argc, char **argv)
{
    if (parse_args(argc, argv) == 1) return depth_usage();

    // if (args.B == NULL && args.tid != -1)
    // return bam2depth(args.idx, args.tid, args.start, args.end, args.strand);

    memory_release();
    
    return 0;
}
