#include "utils.h"

#define compA 0x1
#define compC 0x2
#define compG 0x4
#define compT 0x8
#define compN 0xf
#define compR 0xa  // 10, A or G
#define compY 0x5  // C or T
#define compS 0x6  // G or C
#define compW 0x9  // A or T
#define compK 0x3  // G or T
#define compM 0xc  // 12, A or C
#define compB 0x7  //  C or G or T
#define compD 0xb  // 11, A or G or T
#define compH 0xd  // 13, A or C or T
#define compV 0xe  // 14, A or C or G


uint8_t compDNA_map_table[256] = {
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    1,  2,  4,  8,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  1,  7,  2, 11,  0,  0,  4, 13,  0,  0,  3,  0, 12, 15,  0,
    0,  0, 10,  6,  8,  0, 14,  9,  0,  5,  0,  0,  0,  0,  0,  0,
    0,  1,  7,  2, 11,  0,  0,  4, 13,  0,  0,  3,  0, 12, 15,  0,
    0,  0, 10,  6,  8,  0, 14,  9,  0,  5,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
};

/*
    compDNA_map_table['a'] = 0x1;
    compDNA_map_table['A'] = 0x1;
    compDNA_map_table['c'] = 0x2;
    compDNA_map_table['C'] = 0x2;
    compDNA_map_table['g'] = 0x4;
    compDNA_map_table['G'] = 0x4;
    compDNA_map_table['t'] = 0x8;
    compDNA_map_table['T'] = 0x8;
    compDNA_map_table['n'] = 0xf;
    compDNA_map_table['N'] = 0xf;
    compDNA_map_table['r'] = 0xa;
    compDNA_map_table['R'] = 0xa;
    compDNA_map_table['y'] = 0x5;
    compDNA_map_table['Y'] = 0x5;
    compDNA_map_table['s'] = 0x6;
    compDNA_map_table['S'] = 0x6;
    compDNA_map_table['w'] = 0x9;
    compDNA_map_table['W'] = 0x9;
    compDNA_map_table['k'] = 0x3;
    compDNA_map_table['K'] = 0x3;
    compDNA_map_table['m'] = 0xc;
    compDNA_map_table['M'] = 0xc;
    compDNA_map_table['b'] = 0x7;
    compDNA_map_table['B'] = 0x7;
    compDNA_map_table['d'] = 0xb;
    compDNA_map_table['D'] = 0xb;
    compDNA_map_table['h'] = 0xd;
    compDNA_map_table['H'] = 0xd;
    compDNA_map_table['v'] = 0xe;
    compDNA_map_table['V'] = 0xe;

    int k1, k2;
    for (k1 = 0; k1 < 16; ++k1) {
        for (k2 = 0; k2 < 16; ++k2) {
            fprintf(stderr, "%2d, ", compDNA_map_table[k1*16+k2]);
        }
        fprintf(stderr, "\n");
    }
    exit(1);
*/
const char compDNA_str[] = "\0ACKGYSBTWRDMHVN";

char *compactDNA(const char *a, int l)
{
    int l0 = l/2 + 1 + (l&1);
    char *c = malloc(l0);
    memset(c, 0, l0);
    int i, j = 0;
    for (i = 0; i < l; ++i) {
        uint8_t c1 = compDNA_map_table[(uint8_t)a[i]];
        if (c1 == 0) error("Try to compact unknown base, %c", a[i]);
        
        if (c[j]) {
            c[j] |= (c1 & 0xf);
            j++;
        }
        else {
            c[j] = ((c1&0xf)<<4)&0xf0;
        }
    }
    c[l0] = 0;
    return c;
}

char *compDNA_decode(const char *a)
{
    assert(a);
    int l;
    l = strlen(a);
    int l1 = l*2 + (a[l-1]&0xf);
    char *r = malloc(l1);
    int i;
    int j = 0;
    for (i = 0; i < l; ++i) {
        r[j++] = compDNA_str[((a[i]>>4)&0xf)];
        r[j++] = compDNA_str[(a[i]&0xf)];
    }
    r[j] = '\0';
    //r[l1-1] = '\0';
    return r;
}

int compDNA_hamming_distance(const char *a, const char *b)
{
    assert(a);
    assert(b);

    int l0 = strlen(a);
    int l1 = strlen(b);
    if (l0 != l1) error("Unequal length.");
    int i;
    int e = 0;
    for (i = 0; i < l0; ++i) {
        if (a[i] != b[i]) {
            uint8_t x = a[i] ^ b[i];
            if (x & 0xf ) e++;
            if (x & 0xf0) e++;
        }
    }
      
    return e;
}
