#include "utils.h"
#include "pisa_version.h"
#include <string.h>

typedef struct {
    const char *version;
    const char *date;
    const char *note;
} pisa_release_t;

static const pisa_release_t releases[] = {
    {"v1.8", "2026-06-02", ""},
    {"v1.7", "2026-04-15", ""},
    {"v1.6", "2025-05-29", ""},
    {"v1.5", "2025-03-17", ""},
    {"v1.4", "2025-03-15", ""},
    {"v1.3", "2024-10-10", ""},
    {"v1.2", "2024-09-07", ""},
    {"v1.1", "2024-06-11", ""},
    {"v1.0b", "2023-12-02", "pre-release"},
    {"v1.0a", "2023-10-16", ""},
    {"v0.12", "2022-07-01", ""},
    {"v0.12b", "2022-04-26", ""},
    {"v0.12a", "2022-03-30", ""},
    {"v0.11a", "2022-03-13", ""},
    {"v0.10", "2022-01-06", ""},
    {"v0.10b", "2021-12-09", ""},
    {"v0.10a", "2021-11-08", ""},
    {"v0.9",  "2021-10-14", ""},
    {"v0.8",  "2021-07-20", ""},
    {"v0.7",  "2021-03-26", ""},
    {"v0.6",  "2020-10-29", ""},
    {"v0.5",  "2020-08-27", ""},
    {"v0.4",  "2020-07-13", ""},
    {"v0.3",  "2020-03-26", ""},
    {"0.3",   "2020-03-02", ""},
    {"0.1",   "2019-11-13", ""},
    {"0.0.0.9999", "2019-05-19", "initial"},
    {NULL, NULL, NULL}
};

static const char *cite_bioinf =
    "Shi Q, Liu S, Kristiansen K, Liu L. The FASTQ+ format and PISA. "
    "Bioinformatics. 2022;38(19):4639-4642. doi:10.1093/bioinformatics/btac562";

static const char *cite_yano =
    "Shi Q, Kristiansen K. Spatial dissimilarity analysis in single cell "
    "transcriptomics. Cell Reports Methods. 2025. doi:10.1016/j.crmeth.2025.101141";

static int history_usage()
{
    fprintf(stderr, "# Show PISA release history.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m history\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -c        Show only stable releases (no alpha/beta).\n");
    fprintf(stderr, " -n [N]    Show only the last N releases.\n");
    fprintf(stderr, " -h        Print this help.\n");
    fprintf(stderr, "\n");
    return 1;
}

static int is_stable(const char *v)
{
    /* Treat versions with 'a'/'b' suffix or leading '0.' as pre-release */
    const char *p = v;
    if (*p == '0') return 0;
    for (; *p; p++) {
        if (*p == 'a' || *p == 'b') return 0;
    }
    return 1;
}

int pisa_history(int argc, char *argv[])
{
    int stable_only = 0;
    int last_n = 0;
    int i;

    for (i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0) {
            return history_usage();
        } else if (strcmp(argv[i], "-c") == 0) {
            stable_only = 1;
        } else if (strcmp(argv[i], "-n") == 0 && i + 1 < argc) {
            last_n = atoi(argv[++i]);
        }
    }

    fprintf(stdout, "\n\x1b[1mPISA Release History\x1b[0m\n");
    fprintf(stdout, "Current version: \x1b[36m\x1b[1m%s\x1b[0m\n", PISA_VERSION);
    fprintf(stdout, "Homepage: https://github.com/shiquan/PISA\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "%-16s  %-12s  %s\n", "VERSION", "DATE", "NOTE");
    fprintf(stdout, "----------------  ------------  ----------\n");

    int shown = 0;
    for (i = 0; releases[i].version; i++) {
        if (stable_only && !is_stable(releases[i].version)) continue;
        if (last_n > 0 && shown >= last_n) break;

        const char *marker = "";
        if (strcmp(releases[i].version, PISA_VERSION) == 0)
            marker = "\x1b[32m<< current\x1b[0m";

        fprintf(stdout, "\x1b[1m%-16s\x1b[0m  %-12s  %s %s\n",
                releases[i].version,
                releases[i].date,
                releases[i].note,
                marker);
        shown++;
    }

    fprintf(stdout, "\n");
    fprintf(stdout, "\x1b[1mCitation:\x1b[0m\n");
    fprintf(stdout, "  For PISA v1.4+ (2025--):  %s\n", cite_yano);
    fprintf(stdout, "  For PISA v1.3- (--2024):  %s\n", cite_bioinf);
    fprintf(stdout, "\n");
    fprintf(stdout, "Reference: https://github.com/shiquan/PISA/releases\n");
    fprintf(stdout, "\n");
    return 0;
}
