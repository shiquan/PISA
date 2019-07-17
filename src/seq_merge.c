#include "utils.h"
#include "ksw.h"
#include "htslib/kstring.h"

#define MAX_SCORE_RATIO 0.9f

unsigned char nst_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

static void bwa_fill_scmat(int a, int b, int8_t mat[25])
{
	int i, j, k;
	for (i = k = 0; i < 4; ++i) {
		for (j = 0; j < 4; ++j)
			mat[k++] = i == j? a : -b;
		mat[k++] = -1; // ambiguous base
	}
	for (j = 0; j < 5; ++j) mat[k++] = -1;
}
#define MAX_MIS 3
// not check the edges, **risk**
int check_compl(const char *s1, const char *s2, int l)
{
    int i;
    int m = 0;
    for (i = 0; i < l; ++i) {
        if (nst_nt4_table[(int)s1[i]] + nst_nt4_table[(int)s2[l-i-1]] != 3) m++;
        if (m > MAX_MIS) return 1;
    }
    return 0;
}
uint8_t *nt4_enc(char *s, int l)
{
    uint8_t *e = malloc(l);
    int i;
    for (i = 0; i < l; ++i) e[i] = nst_nt4_table[(int)s[i]];
    return e;
}
int check_overlap(const int lr, const int lq, char const *ref, uint8_t *qry)
{
    uint8_t *s;
    int i, xtra;
    kswr_t r;
    int8_t mat[25];
    bwa_fill_scmat(5, 4, mat);
    s = malloc(lr);
    for (i = 0; i < lr; ++i) {
        int c = ref[i];
        s[i] = c < 0 || c > 127? 4 : c <= 4? c : nst_nt4_table[c];
    }

    xtra = KSW_XSTART | KSW_XSUBO;
    r = ksw_align(lq, qry, lr, s, 5, mat, 2, 17, xtra, 0);
    ++r.qe; ++r.te; // change to the half-close-half-open coordinates
    free(s);   
    if (r.score < 40 || r.qe - r.qb != r.te - r.tb) {
        return -1;
    }

    return r.tb;
}
// adapt from bwa/pemerge.c
// return 1 on unchange, 0 on merged
int merge_paired(const int l_seq1, const int l_seq2, char const *s1, char const *s2, char const *q1, char const *q2, char **_seq, char **_qual)
{
    uint8_t *s[2], *q[2], *seq, *qual;
    int i, xtra, l_seq;
    kswr_t r;
    int8_t mat[25];
    bwa_fill_scmat(5, 4, mat);
    s[0] = malloc(l_seq1); q[0] = malloc(l_seq1);
    s[1] = malloc(l_seq2); q[1] = malloc(l_seq2);
    for (i = 0; i < l_seq1; ++i) {
        int c = s1[i];
        s[0][i] = c < 0 || c > 127? 4 : c <= 4? c : nst_nt4_table[c];
        q[0][i] = q1? q1[i] - 33 : 60;
    }
    for (i = 0; i < l_seq2; ++i) {
        int c = s2[l_seq2 - 1 - i];
        c = c < 0 || c > 127? 4 : c < 4? c : nst_nt4_table[c];
        s[1][i] = c < 4? 3 - c : 4;
        q[1][i] = q2? q2[l_seq2 - 1 - i] - 33 : 60;
    }
    
    xtra = KSW_XSTART | KSW_XSUBO;
    r = ksw_align(l_seq2, s[1], l_seq1, s[0], 5, mat, 2, 17, xtra, 0);
    ++r.qe; ++r.te; // change to the half-close-half-open coordinates

    *_seq = NULL; *_qual = NULL;
    if (r.score < 50) goto pem_ret; // poor alignment
    if (r.tb < r.qb) goto pem_ret;
    if (l_seq1 - r.te > l_seq2 - r.qe) goto pem_ret; // no enough space for the right end
    if ((double)r.score2 / r.score >= MAX_SCORE_RATIO) goto pem_ret; // the second best score is too large
    if (r.qe - r.qb != r.te - r.tb) goto pem_ret; // we do not allow gaps
    
    int l = l_seq1 - (r.tb - r.qb); // length to merge
    l_seq = l_seq1 + l_seq2 - l;
    seq = malloc(l_seq + 1);
    qual = malloc(l_seq + 1);
    memcpy(seq,  s[0], l_seq1); memcpy(seq  + l_seq1, &s[1][l], l_seq2 - l);        
    memcpy(qual, q[0], l_seq1); memcpy(qual + l_seq1, &q[1][l], l_seq2 - l);
    for (i = 0; i < l_seq; ++i) seq[i] = "ACGTN"[(int)seq[i]], qual[i] += 33;
    seq[l_seq] = qual[l_seq] = 0;
    *_seq = (char*)seq; *_qual = (char*)qual;
    free(s[0]); free(s[1]); free(q[0]); free(q[1]);
    return 0;

pem_ret:
    free(s[0]); free(s[1]); free(q[0]); free(q[1]);
    return 1;
}

char *check_circle(char *seq)
{   
    uint8_t *s[2];
    int i, xtra, l_seq;
    kswr_t r;
    int8_t mat[25];
    int seed_length = 20;
    int fragment_limit = 200;
    l_seq = strlen(seq);
    if (l_seq < fragment_limit) return NULL; // too short, only check large fragment
    
    bwa_fill_scmat(5, 4, mat);
    
    s[0] = malloc(seed_length);
    s[1] = malloc(l_seq - seed_length);
    for (i = 0; i < seed_length; ++i) {
        int c = seq[i];
        s[0][i] = c < 0 || c > 127? 4 : c <= 4? c : nst_nt4_table[c];
    }
    for (i = 0; i < l_seq - seed_length; ++i) {
        int c = seq[i+seed_length];
        c = c < 0 || c > 127? 4 : c < 4? c : nst_nt4_table[c];
        s[1][i] = c < 0 || c > 127? 4 : c <= 4? c : nst_nt4_table[c];
    }
    
    xtra = KSW_XSTART | KSW_XSUBO;
    r = ksw_align(seed_length, s[0], l_seq-seed_length, s[1], 5, mat, 2, 17, xtra, 0);
    free(s[0]); free(s[1]);
    if (r.qe - r.qb == r.te - r.tb && r.qe - r.qb == seed_length-1 && r.score >= 40) {
        /*
        int j, k;
        int mis = 0;        
        for (j = r.te, k = r.qe; j < l_seq; ++j, ++k) 
            if (s[j] != s[k]) mis++;
        if ((float)mis/(l_seq-j) > 0.1) return NULL;
        */
        kstring_t str = {0,0,0};
        kputsn(seq, r.te+1, &str);
        kputs("", &str);
        return str.s;
    }
    return NULL;
}

/*
int main(int argc, char **argv)
{
    if (argc != 2) error("check sequence");
    char *seq = strdup(argv[1]);
    printf("%s\n", seq);
    char *s = check_circle(seq);
    printf("%s\n", s);
    free(seq); free(s);
    return 0;
}
*/
