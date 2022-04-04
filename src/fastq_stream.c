#include "utils.h"
#include "fastq.h"
#include "number.h"
#include "dict.h"
#include "fastq.h"
#include "read_tags.h"
#include "htslib/thread_pool.h"
/*
#ifdef _OPENMP
#include <omp.h>
#endif
*/
static struct args {
    const char *input_fname;
    const char *output_fname;
    const char *report_fname;
    const char *script;
    const char *tags_str;
    const char *tempdir;
    struct dict *tags;
    
    int min_reads_per_block;
    int max_reads_per_block;
    int keep_processed;
    struct fastq_handler *fastq;
    int n_thread;

    char *run_script;
    int no_warnings;
    
    // In default, temp files will be deleted after process, but enable this flag will keep them.
    // However, program will create one temp dir for each read block which means it will be
    // harmful to the file system if much blocks are processed. For example, to stream over 10K
    // cells will create over 10K folders in the temp directory. Consider this, PISA will check
    // the number of fastq blocks when you enable this flag. If over 1000 blocks are processed,
    // the keep_temp will be disabled.
    
    int keep_temp;
    int stream_input_fasta;
    FILE *fout;    
} args = {
    .input_fname  = NULL,
    .output_fname = NULL,
    .report_fname = NULL,
    .script       = NULL,
    .tags_str     = NULL,
    .tempdir      = "_PISA_stream_tempdir",
    .tags         = NULL,
    .min_reads_per_block = 2,
    .max_reads_per_block = 8000,
    .keep_processed = 0,
    .fastq        = NULL,
    .n_thread     = 1,
    .run_script   = NULL,
    .keep_temp    = 0,
    .stream_input_fasta = 0,
    .fout         = NULL,
};

static void memory_release()
{
    if (args.run_script) free(args.run_script);
    fastq_handler_destory(args.fastq);    
}

static int parse_args(int argc, char **argv)
{
    const char *thread = NULL;
    const char *min = NULL;
    const char *max = NULL;
    int i;

    for (i = 1; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;

        if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return 1;
        if (strcmp(a, "-o") == 0) var = &args.output_fname;
        else if (strcmp(a, "-t") == 0) var = &thread;
        else if (strcmp(a, "-tags") == 0) var = &args.tags_str;
        else if (strcmp(a, "-tmpdir") == 0) var = &args.tempdir;
        else if (strcmp(a, "-script") == 0) var = &args.script;
        else if (strcmp(a, "-min") == 0) var = &min;
        else if (strcmp(a, "-max") == 0) var = &max;
        else if (strcmp(a, "-nw") == 0) {
            args.no_warnings = 1;
            continue;
        }
        else if (strcmp(a, "-keep-tmp") == 0) {
            args.keep_temp = 1;
            continue;
        }
        else if (strcmp(a, "-fa") == 0) {
            args.stream_input_fasta = 1;
            continue;
        }
        else if (strcmp(a, "-keep") == 0) {
            args.keep_processed = 1;
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
    if (args.script == NULL) error("No -script specified.");
    args.fout = args.output_fname == NULL ? stdout : fopen(args.output_fname, "w");
    if (args.fout == NULL) error("%s : %s.", args.output_fname, strerror(errno));
    
    if (args.tags_str == NULL) error("-tag must be set.");
    args.tags = str2tag(args.tags_str);
    if (thread) args.n_thread = str2int(thread);
    // if (file_thread) args.n_thread = str2int(file_thread);
    if (args.n_thread < 1) args.n_thread = 1;

    if (min) args.min_reads_per_block = str2int(min);
    if (max) args.max_reads_per_block = str2int(max);
    
    struct stat st = {0};

    if (stat(args.tempdir, &st) == -1) {
        mkdir(args.tempdir, 0700);
    }

    args.fastq = fastq_handler_init(args.input_fname, NULL, 0, 0);
    if (args.fastq == NULL) error("%s : %s.", args.input_fname, strerror(errno));
    return 0;    
}

static struct bseq_pool *stream_process(const char *run_script, struct bseq_pool *in, char *tmpdir, char *unique_block_name)
{
    // char cwd[PATH_MAX];
    // getcwd(cwd, sizeof(cwd));
    struct stat st = {0};
    if (stat(tmpdir, &st) == -1) mkdir(tmpdir, 0700);

    // chdir(tmpdir);

    kstring_t script = {0,0,0};
    kstring_t input = {0,0,0};

    kputs(tmpdir, &input);
    kputs("/_block.fq", &input);
    
    // create tmp fastq input, setup $FQ envion parameter(s)
    bseq_pool_write_file(in, input.s);
    free(input.s);
    kputs("cd ", &script);
    kputs(tmpdir, &script);
    kputs("; ", &script);
    kputs("export FQ=\"./_block.fq\"; ", &script);
    ksprintf(&script, "export UBI=\"%s\"; ", unique_block_name);
    
    kputs(run_script, &script);

    FILE *fp;
    fp = popen(script.s, "r");
    if (fp == NULL) {
        if (args.no_warnings) return NULL;

        warnings("Failed to process read block %s : %s.", unique_block_name, strerror(errno));
        // out of tmp dir
        return NULL;
    }
    struct bseq_pool *out = bseq_pool_cache_fp(fp, 0);
    pclose(fp);

    free(script.s);
    
    // back to workdir
    // chdir(cwd); 
    return out;
}

// Run script example.
// "seqtk seq -A $FQ > input.fa && cap3 input.fa -r 0 -a 2 && seqtk rename $UBI_ input.fa.cap.contigs 1> contigs.fa && cat contigs.fa test.fa.cap.singlets"

// struct bseq_pool *stream_script_stdout(const char *script, struct bseq_pool *in)
char *stream_script_format(const char *script)
{
    kstring_t str = {0,0,0};
    
    struct stat sb;
    int file_check = stat(script, &sb);
    if (file_check == 0) { // a script file
        char actualpath [PATH_MAX+1];
        char *path = realpath(script, actualpath);
        if (path == NULL) error("Unable to define the path of script, %s", script);
        ksprintf(&str, "sh %s", actualpath);
        LOG_print("script path: %s", actualpath);
    } else {
        // code
        kputs(script, &str);
    }

    return str.s;
}
extern int fastq_stream_usage();

static void write_out(struct bseq_pool *p)
{
    if (p == NULL) return;        
    bseq_pool_write_fp(p, args.fout);

    if (p->opts) {
        int i;
        char **vals = p->opts;
        int n = dict_size(args.tags);
        for (i = 0; i < n; ++i)
            if (vals[i]) free(vals[i]);
        free(vals);
    }
    bseq_pool_destroy(p);
}
static void *run_it(void *_data)
{
    struct bseq_pool *p = (struct bseq_pool*)_data;
    
    if (p->n == 1 || p->n < args.min_reads_per_block) {
        if (args.keep_processed == 0) return NULL;
        if (args.stream_input_fasta == 1) p->force_fasta = 1;
        return p;
    }

    // make unique block barcode [0-9A-Za-z_]
    char *ubi = vals2str(p->opts, dict_size(args.tags));

    kstring_t tempdir0 = {0,0,0};
    kputs(args.tempdir, &tempdir0);
    if (tempdir0.s[tempdir0.l-1] != '/') kputc('/', &tempdir0);
    kputs(ubi, &tempdir0);
    struct bseq_pool *ret_p = stream_process(args.run_script, p, tempdir0.s, ubi);
    
    if (ret_p == NULL) {
        if (args.no_warnings == 0) warnings("Block %s has zero output.", ubi);
    } else {
        if (args.stream_input_fasta == 1) ret_p->force_fasta = 1;
        // update names
        char **vals = p->opts;
        int i;
        for (i = 0; i < ret_p->n; ++i) {
            struct bseq *b = &ret_p->s[i];
            char *new = fname_update_tags(b->n0.s, args.tags, vals);
            if (new) {
                b->n0.l = 0;
                kputs(new, &b->n0);
                free(new);
            }
        }
    }

    int i;
    char **vals = p->opts;
    int n = dict_size(args.tags);
    for (i = 0; i < n; ++i)
        if (vals[i]) free(vals[i]);
    free(vals);
    bseq_pool_destroy(p);

    free(ubi);

    if (args.keep_temp == 0) {
        kstring_t str = {0,0,0};
        kputs("rm -rf ", &str);
        kputs(tempdir0.s, &str);
        
        if (system(str.s) == -1) warnings("Failed to run %s", str.s);
        free(str.s);
    }
    free(tempdir0.s);

    return ret_p;
}

int fastq_stream(int argc, char **argv)
{
    double t_real;
    t_real = realtime();
    
    if (parse_args(argc, argv)) return fastq_stream_usage();
    
    args.run_script = stream_script_format(args.script);
    if (args.run_script == NULL) error("Empty run script?");
    
    int n_block = 0;

    hts_tpool *p = hts_tpool_init(args.n_thread);
    hts_tpool_process *q = hts_tpool_process_init(p, args.n_thread*2, 0);
    hts_tpool_result *r;
    
    for (;;) {
        if (args.fastq->closed == 1) break;
        struct bseq_pool *b0;
        b0 = fastq_read_block(args.fastq, args.tags);
        if (b0 == NULL) break;
        n_block++;
        
        if (n_block > 1000 && args.keep_temp ==1) {
            warnings("Too much temp files created. Auto disable -keep-tmp.");
            args.keep_temp = 0;   
        }
        
        int block;
        do {
            block = hts_tpool_dispatch2(p, q, run_it, b0, 1);
            if ((r = hts_tpool_next_result(q))) {
                struct bseq_pool *d = (struct bseq_pool*)hts_tpool_result_data(r);
                write_out(d);
                hts_tpool_delete_result(r, 0);
            }
        }
        while (block == -1);
    }
    
    hts_tpool_process_flush(q);
    
    while ((r = hts_tpool_next_result(q))) {
        struct bseq_pool *d = (struct bseq_pool*)hts_tpool_result_data(r);
        write_out(d);
        hts_tpool_delete_result(r, 0);
    }
    hts_tpool_process_destroy(q);
    hts_tpool_destroy(p);

    memory_release();

    LOG_print("Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB.", realtime() - t_real, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
    return 0;
}


/*
int fastq_stream0(int argc, char **argv)
{
    double t_real;
    t_real = realtime();
    
    if (parse_args(argc, argv)) return fastq_stream_usage();
    
    char *run_script = stream_script_format(args.script);
    if (run_script == NULL) error("Empty run script?");
    
    int n_block = 0;
    int n = dict_size(args.tags);

    
#pragma omp parallel shared(n_block) num_threads(args.n_thread)
    for (;;) {
        if (args.fastq->closed == 1) break;
        struct bseq_pool *p;
#pragma omp critical (read)
        p = fastq_read_block(args.fastq, args.tags);
        if (p == NULL) break;
#pragma omp atomic
        n_block++;
        
        if (n_block > 1000 && args.keep_temp ==1) {
            warnings("Too much temp files created. Auto disable -keep-tmp.");
            args.keep_temp = 0;   
        }
        
        if (p->n == 1 || p->n < args.min_reads_per_block) {
            if (args.keep_processed == 0) continue;
            if (args.stream_input_fasta == 1) p->force_fasta = 1;

#pragma omp critical (write)
            bseq_pool_write_fp(p, args.fout);
            
            int i;
            char **vals = p->opts;
            for (i = 0; i < n; ++i)
                if (vals[i]) free(vals[i]);
            free(vals);
                     
            bseq_pool_destroy(p);
            continue;
        }

        // make unique block barcode [0-9A-Za-z_]
        char *ubi = vals2str(p->opts, n);

        kstring_t tempdir0 = {0,0,0};
        kputs(args.tempdir, &tempdir0);
        if (tempdir0.s[tempdir0.l-1] != '/') kputc('/', &tempdir0);
        kputs(ubi, &tempdir0);
        struct bseq_pool *ret_p = stream_process(run_script, p, tempdir0.s, ubi);
        
        if (ret_p == NULL) {
            if (args.no_warnings == 0) warnings("Block %s has zero output.", ubi);
        } else {
            if (args.stream_input_fasta == 1) ret_p->force_fasta = 1;
            // update names
            char **vals = p->opts;
            int i;
            for (i = 0; i < ret_p->n; ++i) {
                struct bseq *b = &ret_p->s[i];
                char *new = fname_update_tags(b->n0.s, args.tags, vals);
                if (new) {
                    b->n0.l = 0;
                    kputs(new, &b->n0);
                    free(new);
                }
            }
#pragma omp critical (write)
            bseq_pool_write_fp(ret_p, args.fout);
            bseq_pool_destroy(ret_p);
        }
        
        free(ubi);
        
        char **vals = p->opts;
        int i;
        for (i = 0; i < n; ++i) free(vals[i]);
        free(vals);
        bseq_pool_destroy(p);

        if (args.keep_temp == 0) {
            kstring_t str = {0,0,0};
            kputs("rm -rf ", &str);
            kputs(tempdir0.s, &str);
            
            system(str.s);
            free(str.s);
        }
        free(tempdir0.s);
    }

    free(run_script);
    memory_release();

    LOG_print("Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB.", realtime() - t_real, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
    return 0;
}
*/
