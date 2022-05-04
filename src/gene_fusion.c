// gene_fusion.c -- intepret gene fusions based on UMI tagged cDNA library.
//
//
#include "utils.h"
#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "number.h"
#include "dict.h"

static struct args {
    const char *input_fname;
    int qual_thres;
    const char *gene_tag;
    const char *cell_tag;
    const char *umi_tag;
    const char *cell_list;
    const char *output_fname;
    const char *gtf_fname;
    const char *fs_tag;
    struct dict *cell_barcodes;
    struct dict *genes;
    int file_th;
    int mmgu;
} args = {
    .input_fname   = NULL,
    .qual_thres    = 0,
    .gene_tag      = "GN",
    .cell_tag      = "CB",
    .umi_tag       = "UB",
    .fs_tag        = "FS",
    .cell_list     = NULL,
    .output_fname  = NULL,
    .gtf_fname     = NULL,
    .genes         = NULL,
    .cell_barcodes = NULL,
    .file_th       = 4,
    .mmgu          = 3,
};

static int fixed_barcodes = 0;

static int fusion_only = 0;
static int empty_barcodes = 1;

struct barcode_gene {
    int n; // genes, usually be 1
    int *gene;
};

static void memory_release()
{
    int i;
    for (i = 0; i < dict_size(args.cell_barcodes); ++i) {
        struct dict *umi = dict_query_value(args.cell_barcodes, i);
        if (!umi) continue;
        int j;
        for (j = 0; j < dict_size(umi); ++j) {
            struct barcode_gene * bg = dict_query_value(umi, j);
            if (!bg) continue;
            free(bg->gene);
            free(bg);
        }
        dict_destroy(umi);
    }
    dict_destroy(args.cell_barcodes);
    dict_destroy(args.genes);
}
int build_barcode_gene_table()
{
    htsFile *fp = hts_open(args.input_fname, "r");
    bam_hdr_t *hdr = sam_hdr_read(fp);
    
    hts_set_threads(fp, args.file_th);

    kstring_t str = {0,0,0};
    bam1_t *b = bam_init1();
    int ret;
    for (;;) {
        ret = sam_read1(fp, hdr, b);
        if (ret < 0) break;
        
        if (b->core.flag & BAM_FQCFAIL) continue;
        if (b->core.flag & BAM_FSECONDARY) continue;
        if (b->core.qual < args.qual_thres) continue;
        
        uint8_t *data = bam_aux_get(b, args.gene_tag);
        if (!data) continue; // if no gene tag

        empty_barcodes = 0;
        
        str.l = 0;
        kputs((char*)(data+1), &str);
        int i0 = 0;
        int overlap_gene = 0; // if gene overlaped, skip them
        for (i0 = 0; i0 < str.l; ++i0) {
            if (str.s[i0] == ';') {
                overlap_gene = 1;
                break;
            }
        }

        if (overlap_gene == 1) continue;
        
        int gene_idx = dict_push(args.genes, (char*)(data+1));
        
        data = bam_aux_get(b, args.cell_tag);
        if (!data) continue;
        
        int cell_idx = dict_query(args.cell_barcodes, (char*)(data+1));
        if (fixed_barcodes == 1 && cell_idx == -1) continue;
        
        cell_idx = dict_push(args.cell_barcodes, (char*)(data+1));
        struct dict *umi = dict_query_value(args.cell_barcodes, cell_idx);

        if (umi == NULL) {
            umi = dict_init();
            dict_set_value(umi);
            dict_assign_value(args.cell_barcodes, cell_idx, umi);
        }

        data = bam_aux_get(b, args.umi_tag);
        if (!data) continue;
        
        int umi_idx = dict_push(umi, (char*)(data+1));
        
        struct barcode_gene *bg = dict_query_value(umi, umi_idx);
        if (bg == NULL) {
            bg = malloc(sizeof(struct barcode_gene));
            bg->n = 0;
            bg->gene = NULL;
            dict_assign_value(umi, umi_idx, bg);
        }

        // In case duplicate gene name for a umi
        int i;
        for (i = 0; i < bg->n; ++i) {
            if (gene_idx == bg->gene[i]) break;
        }

        if (i == bg->n) { // not found
            if (bg->n == 0) {
                bg->gene = malloc(sizeof(int));
                bg->gene[0] = gene_idx;
                bg->n = 1;
            }
            else {
                bg->gene = realloc(bg->gene, sizeof(int)*(bg->n+1));
                bg->gene[bg->n] = gene_idx;
                bg->n++;
            }
        }
    }

    if (str.m) free(str.s);
    sam_close(fp);
    bam_hdr_destroy(hdr);
    bam_destroy1(b);

    return empty_barcodes;
}

// Remove non fusion records, for query effcient.
int slim_barcode_gene_table()
{
    int i;
    int fusions = 0;
    for (i = 0; i < dict_size(args.cell_barcodes); ++i) {

        int delete_this_cell = 1;
        struct dict *umi = dict_query_value(args.cell_barcodes, i);
        int j;
        for (j = 0; j < dict_size(umi); ++j) {
            struct barcode_gene *bg = dict_query_value(umi, j);
            assert(bg); // at this memonet, all value should be assigned

            if (bg->n == 1) {
                free(bg->gene);
                free(bg);
                dict_assign_value(umi, j, NULL);
            }
            
            // too much hits for a umi, perhaps low quality, skip it then
            
            else if (bg->n > args.mmgu) {
                free(bg->gene);
                free(bg);
                dict_assign_value(umi, j, NULL);
            }

            else if (bg->n > 1) {
                delete_this_cell = 0; // disable flag
                fusions++;
            }
        }

        if (delete_this_cell) {
            dict_destroy(umi);
            dict_assign_value(args.cell_barcodes,i, NULL);
        }
    }

    return fusions;
}

static int parse_args(int argc, const char **argv)
{
    int i;
    const char *file_thread = NULL;
    const char *map_qual = NULL;
    const char *max = NULL;
    
    for (i = 1; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;

        // common options
        if (strcmp(a, "-o") == 0 ) var = &args.output_fname;
        else if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return 1;
        else if (strcmp(a, "-@") == 0) var = &file_thread;
        else if (strcmp(a, "-q") == 0) var = &map_qual;
        else if (strcmp(a, "-gn") == 0) var = &args.gene_tag;
        else if (strcmp(a, "-cb") == 0) var = &args.cell_tag;
        else if (strcmp(a, "-umi") == 0) var = &args.umi_tag;
        else if (strcmp(a, "-list") == 0) var = &args.cell_list;
        else if (strcmp(a, "-m") == 0) var = &max;
        else if (strcmp(a, "-fusion-only") == 0) {
            fusion_only = 1;
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


    if (args.input_fname == NULL) error("Input bam must be set.");

    if (file_thread) args.file_th = str2int(file_thread);

    if (max) args.mmgu = str2int(max);

    if (map_qual) args.qual_thres = str2int(map_qual);

    args.cell_barcodes = dict_init();
    args.genes = dict_init();
    
    dict_set_value(args.cell_barcodes);

    if (args.cell_list) {
        dict_read(args.cell_barcodes, args.cell_list, 0);
        fixed_barcodes = 1;
    }
    
    // Check if bam indexed, if not index it.

    // Build cell barcode-UMI and mapped gene table.
    LOG_print("Building barcode table..");
    if (build_barcode_gene_table() == 1)
        error("No found gene tag in the bam. Did you `PISA anno` your bam file?");

    return 0;
}
extern int gene_fusion_usage();

int gene_fusion(int argc, const char **argv)
{
    double t_real;
    t_real = realtime();
    
    if (parse_args(argc, argv) == 1) return gene_fusion_usage();
    LOG_print("Slimming barcode table..");
    if (slim_barcode_gene_table() == 0) {
        LOG_print("0 fusion reads detected.");
        LOG_print("Real time: %.3f sec; CPU: %.3f sec.", realtime() - t_real, cputime());

        return 0;
    }
    LOG_print("Start annotation..");
    
    // TODO: breakpoint finding.
    
    // Annotate fusion record into bam file.
    htsFile *fp = hts_open(args.input_fname, "r");
    bam_hdr_t *hdr = sam_hdr_read(fp);
    htsFile *out = hts_open(args.output_fname, "wb");

    if (sam_hdr_write(out, hdr)) error("Failed to write SAM header.");

    hts_set_threads(fp, args.file_th);

    kstring_t str = {0,0,0};
    uint64_t fusion_reads = 0;
    bam1_t *b = bam_init1();
    int ret;
    int fusion;

    for (;;) {
        ret = sam_read1(fp, hdr, b);
        if (ret < 0) break; 

        fusion = 0; // reset fusion flag
        
        if (b->core.flag & BAM_FQCFAIL) goto write_file;
        if (b->core.flag & BAM_FSECONDARY) goto write_file;
        if (b->core.qual < args.qual_thres) goto write_file;

        uint8_t *data = bam_aux_get(b, args.gene_tag);
        if (!data) goto write_file;
        
        data = bam_aux_get(b, args.cell_tag);
        if (!data) goto write_file;
        
        int cell_idx = dict_query(args.cell_barcodes, (char*)(data+1));
        if (cell_idx == -1) goto write_file;
        
        struct dict *umi = dict_query_value(args.cell_barcodes, cell_idx);
        
        // Since table has been slimmed, non-fusion records have be delected.
        if (umi == NULL) goto write_file;
        data = bam_aux_get(b, args.umi_tag);
        if (!data) goto write_file;
        
        int umi_idx = dict_query(umi, (char*)(data+1));

        if (umi_idx == -1) goto write_file;
        
        struct barcode_gene *bg = dict_query_value(umi, umi_idx);
        if (bg == NULL) goto write_file;
        
        assert(bg->n > 1); // Map more than 1 gene.

        str.l = 0;
        int i, j;
        for (i = 0; i < bg->n; ++i) {
            for (j = i+1; j < bg->n; ++j) {
                if (str.l) kputc(';', &str);
                kputs(dict_name(args.genes,bg->gene[i]), &str);
                kputc('|', &str);
                kputs(dict_name(args.genes,bg->gene[j]), &str);
                kputs("", &str);
            }
        }
        
        if (str.l) {
            bam_aux_append(b, args.fs_tag, 'Z', str.l+1, (uint8_t*)str.s);
            fusion_reads++;

            fusion = 1;   // set flag
        }
        
        
    write_file:
        if (fusion_only == 1 && fusion == 0) continue;
        if (sam_write1(out, hdr, b) == -1) error("Failed to write SAM."); 
    }
    
    if (str.m) free(str.s);
    bam_destroy1(b);
    bam_hdr_destroy(hdr);
    sam_close(fp);
    sam_close(out);
    
    memory_release();
    
    LOG_print("%ld fusion reads detected.", fusion_reads);
    LOG_print("Real time: %.3f sec; CPU: %.3f sec.", realtime() - t_real, cputime());
    return 0;
}
