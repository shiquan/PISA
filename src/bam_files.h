#ifndef BAM_FILES_H
#define BAM_FILES_H

#include "utils.h"
#include "htslib/hts.h"
#include "htslib/sam.h"

enum bam_file_state {
    file_not_open, // in case a lot of files need to open in parallel, may exceed the limitation
    file_is_open,
    file_closed,
};

struct bam_file {
    char *fname;
    htsFile *fp;
    bam_hdr_t *hdr;
    char *alias;
    enum bam_file_state state;
};

struct bam_files {
    int n;
    int i; // current file
    int n_thread;
    struct bam_file *files;
};

struct bam_files *init_bam_line(const char *bams, int n_thread);
struct bam_files *init_bam_list(const char *file_list, int n_thread);
void close_bam_files(struct bam_files *files);

int read_bam_files(struct bam_files *files, bam1_t *b);

bam_hdr_t *get_hdr(struct bam_files *files);
char *get_alias(struct bam_files *files);
char *get_fname(struct bam_files *files);

#endif
