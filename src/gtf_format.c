#include "utils.h"
#include "gtf.h"
#include "htslib/kstring.h"
#include "bed.h"

int gtf_format_usage()
{
    fprintf(stderr, "# Format and reorder GTF file.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m gtffmt in.gtf\n");
    fprintf(stderr, "\nOptions:\n");
    fprintf(stderr, " -o      [FILE]    Output GTF file\n");
    fprintf(stderr, " -f                Only export gene, transcript, exon and CDS records.\n");
    fprintf(stderr, " -key    [all]     Export selected keys in optional fields.\n");
    fprintf(stderr, " -report [stderr]  Summary report file.");
    fprintf(stderr, "\n");
    return 1;
}

int gtf_format(int argc, char **argv)
{
    if (argc == 1) return gtf_format_usage();
    int only_exon = 0;
    int i;
    const char *input_fname = NULL;
    const char *output_fname = NULL;
    const char *key_str = NULL;
    const char *report_fname = NULL;
    
    for (i = 1; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;

        if (strcmp(a, "-h") == 0) return gtf_format_usage();
        
        if (strcmp(a, "-o") == 0 || strcmp(a, "-out") == 0)
            var = &output_fname;
        else if (strcmp(a, "-key") == 0)
            var = &key_str;
        else if (strcmp(a, "-f") == 0) {
            only_exon = 1;
            continue;
        }
        else if (strcmp(a, "-report") == 0)
            var = &report_fname;
        
        if (var != 0) {
            if (i == argc) error("Miss an argument after %s.", a);
            *var = argv[i++];
            continue;
        }
        
        if (input_fname == NULL) {
            input_fname = a;
            continue;
        }
        error("Unknown argument, %s", a);
    }

    // if (output_fname == NULL) error("Please specify output file by -o.");
    if (input_fname  == NULL) error("No input specified.");

    struct dict *keys = NULL;
    if (key_str) {
        int *s;
        int n;
        kstring_t str = {0,0,0};
        kputs(key_str, &str);
        s = ksplit(&str,',', &n);
        int i;
        keys = dict_init();
        for (i = 0; i < n; ++i) {
            char *key = str.s + s[i];
            dict_push(keys, key);
        }
    }
    struct gtf_spec *G = gtf_read(input_fname, only_exon);
    gtf_dump(G, output_fname, keys);

    if (output_fname) {
        FILE *fr = report_fname == NULL ? stderr : fopen(report_fname, "w");
        if (fr == NULL) error("%s : %s.", report_fname, strerror(errno));

        struct bed_spec *genes = bed_spec_init();
        struct bed_spec *exons = bed_spec_init();
        
        int n_gene_0 = 0;
        int n_gene_1 = 0;
        int n_coding = 0;
        int n_trans = 0;
        uint32_t gene_cov_0 = 0;
        uint32_t gene_cov_1 = 0;
        uint32_t exon_cov = 0;
        int i;
        for (i = 0; i < dict_size(G->name); ++i) {
            struct gtf_ctg *ctg = dict_query_value(G->name,i);
            int j;
            for (j = 0; j < ctg->n_gtf; ++j) {
                struct gtf *gtf = ctg->gtf[j];
                if (gtf->coding) n_coding++;

                if (gtf->strand == 0) n_gene_0++;
                else n_gene_1++;
                
                // push gene
                bed_spec_push0(genes, GTF_seqname(G, gtf->seqname), gtf->start, gtf->end, gtf->strand, NULL, NULL);
                
                int k;
                for (k = 0; k < gtf->n_gtf; ++k) { // trans
                    struct gtf *trans = gtf->gtf[k];
                    n_trans += trans->n_gtf;
                    
                    int l;
                    for (l = 0; l < trans->n_gtf; ++l) { // exon
                        struct gtf *exon = trans->gtf[l];
                        // push exon
                        bed_spec_push0(exons, GTF_seqname(G, exon->seqname), exon->start, exon->end, exon->strand, NULL, NULL);
                    }
                }
            }
        }

        bed_spec_merge0(genes, 1, 0);
        for (i = 0; i < genes->n; ++i) {
            struct bed *b = &genes->bed[i];
            if (b->strand == 1) gene_cov_1 += (b->end - b->start);
            else gene_cov_0 += (b->end - b->start);
        }

        bed_spec_merge0(exons, 1, 0);
        for (i = 0; i < exons->n; ++i) {
            struct bed *b = &exons->bed[i];
            exon_cov += (b->end - b->start);
        }

        bed_spec_destroy(genes);
        bed_spec_destroy(exons);
        
        fprintf(fr, "Genes (forward strand),%d\n", n_gene_0);
        fprintf(fr, "Genes (reverse strand),%d\n", n_gene_1);
        fprintf(fr, "Coding genes,%d\n", n_coding);
        fprintf(fr, "Transcripts,%d\n", n_trans);
        fprintf(fr, "Gene coverage (forward),%d\n", gene_cov_0);
        fprintf(fr, "Gene coverage (reverse),%d\n", gene_cov_1);
        fprintf(fr, "Exon coverage,%d\n", exon_cov);
        uint32_t intron_cov = gene_cov_0 + gene_cov_1 - exon_cov;
        fprintf(fr, "Intron coverage,%d (%.2f%%)\n", intron_cov, (float)intron_cov/(gene_cov_0+gene_cov_1)*100);

        fclose(fr);
    }
    
    gtf_destroy(G);
    if (keys) dict_destroy(keys);
    return 0;
}
