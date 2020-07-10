#include "utils.h"

#define compA 0x1
#define compC 0x2
#define compG 0x4
#define compT 0x8
#define compN 0xf

const uint8_t map_compDNA_table[256] = {
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    1, 2, 4, 8, 15,15,15,15, 15,15,15,15, 15, 0 /*=*/,15,15,
    15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
    15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
    15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
    15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,

    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
};

const char compDNA_str[] = "=ACMGRSVTWYHKDBN";

const int compDNA_int[] = { 4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4 };


char *compactDNA(const char *a, int l)
{
    int l0 = l/2 + l&1;
    char *c = malloc(l0);
    memset(c, 0, l0);
    int i, j = 0;
    for (i = 0; i < l; ++i) {
        uint8_t c1 = map_compDNA_table[a[i]];
        if (c[j]) {    
            c[j] |= (c1 & 0xf);
            j++;
        }
        else {
            c[j] = c1 & 0xf;
            c[j] <<= 4;
        }
    }
    return c;
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
            uint8_t x = a[i] & b[i];
            if (x & 0xf) e++;
            if (x>>4) e++;
        }
    }
    return e;
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
        r[j++] = compDNA_str[a[i] >>4];
        r[j++] = compDNA_str[a[i]&0xf];
    }
    r[l1-1] = '\0';
    return r;
}
