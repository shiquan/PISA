#include "htslib/kstring.h"
#include "bam_files.h"
#include "htslib/sam.h"
#include "htslib/hts.h"

struct bam_files *init_bam_line(const char *bams, int n_thread)
{
    struct bam_files *files = malloc(sizeof(struct bam_files));
    memset(files, 0, sizeof(struct bam_files));

    files->n_thread = n_thread;
    
    kstring_t str = {0,0,0};
    kputs(bams, &str);
    int n;
    int *s = ksplit(&str, ',', &n);
    
    if (n == 1) {
        files->n = 1;
        files->files = malloc(sizeof(struct bam_file) *1);
        struct bam_file *file = &files->files[0];
        memset(file, 0, sizeof(struct bam_file));
        file->fname = strdup(bams);
        file->fp = hts_open(file->fname, "r");
        if (file->fp == NULL) error("%s : %s.", file->fname, strerror(errno));
        file->state = file_is_open;
        htsFormat type = *hts_get_format(file->fp);
        if (type.format != bam && type.format != sam)
            error("Unsupported input format, only support BAM/SAM/CRAM format.");

        if (files->n_thread>1)
            hts_set_threads(file->fp, files->n_thread);
        
        file->hdr = sam_hdr_read(file->fp);
        if (file->hdr == NULL) error("Failed to open bam header of %s", file->fname);
    }
    else {
        files->n = n;
        files->files = malloc(sizeof(struct bam_file) * n);
        
        int i;
        for (i = 0; i < n; ++i) {
            struct bam_file *file = &files->files[i];
            memset(file, 0, sizeof(struct bam_file));
            file->state = file_not_open;
            file->fname = strdup(str.s+s[i]);
        }
    }

    free(str.s);
    free(s);
    
    return files;
}

struct bam_files *init_bam_list(const char *file_list, int n_thread)
{
    int n;
    char **list = hts_readlist(file_list, 1, &n);
    if (n == 0) error("Empty list. %s", file_list);

    struct bam_files *files = malloc(sizeof(struct bam_files));
    memset(files, 0, sizeof(struct bam_files));
    files->n_thread = n_thread;
    files->n = n;
    files->files = malloc(sizeof(struct bam_file)*n);
    
    kstring_t str = {0,0,0};
    
    int i;
    for (i = 0; i < n; ++i) {
        str.l = 0;
        kputs(list[i], &str);
        free(list[i]);
        
        if (str.l == 0) continue;
        if (str.s[0] == '#') continue;
        
        int c = 0;
        int *s = ksplit(&str, '\t', &c);
        if (c == 0) continue; // empty list
        struct bam_file *file = &files->files[i];
        memset(file, 0, sizeof(struct bam_file));

        file->state = file_not_open;
        file->fname = strdup(str.s);
        
        if (c > 1) file->alias = strdup(str.s+s[1]);

        free(s);
    }
    free(list);
    if (str.m) free(str.s);
    return files;
}

void close_bam_files(struct bam_files *files)
{
    int i;
    for (i = 0; i < files->n; ++i) {
        struct bam_file *file = &files->files[i];
        if (file->state == file_is_open) {
            bam_hdr_destroy(file->hdr);
            sam_close(file->fp);
        }
        free(file->fname);
        if (file->alias) free(file->alias);
    }
    free(files->files);
    free(files);
}

int read_bam_files(struct bam_files *files, bam1_t *b)
{
    if (files->i == files->n) return -1; // for muti-threads
    
    struct bam_file *file = &files->files[files->i];

    while (file->state == file_closed && files->i < files->n) {
        files->i++;
        if (files->i == files->n) return -1;
    }
    
    if (file->state == file_not_open) {
        file->fp = hts_open(file->fname, "r");
        if (file->fp == NULL) error("%s : %s.", file->fname, strerror(errno));
        file->state = file_is_open;
        htsFormat type = *hts_get_format(file->fp);
        if (type.format != bam && type.format != sam)
            error("Unsupported input format, only support BAM/SAM/CRAM format.");

        if (files->n_thread > 1)
            hts_set_threads(file->fp, files->n_thread);

        file->hdr = sam_hdr_read(file->fp);
        if (file->hdr == NULL) error("Failed to open bam header of %s", file->fname);
    }
    
    int ret;
    ret = sam_read1(file->fp, file->hdr, b);
    if (ret < 0) {
        if (file->state == file_is_open) {
            bam_hdr_destroy(file->hdr);
            sam_close(file->fp);
            file->state = file_closed;
        }
        
        files->i++;
        if (files->i == files->n) return ret;
        return read_bam_files(files, b);
    }
    
    return ret;
}

bam_hdr_t *get_hdr(struct bam_files *files)
{
    return files->files[files->i].hdr;
}

char *get_alias(struct bam_files *files)
{
    return files->files[files->i].alias;
}
char *get_fname(struct bam_files *files)
{
    return files->files[files->i].fname;
}
