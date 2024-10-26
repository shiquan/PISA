#include "utils.h"
#include "read_anno.h"

static enum exon_type RE_map[256] = {
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
    0,  7,  8,  3,  0,  2,  0,  0,  0, 10,  0,  0,  0,  0,  5,  0, 
    0,  0,  4,  1,  0,  0,  9,  0,  6,  0,  0,  0,  0,  0,  0,  0, 
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
};

enum exon_type RE_type_map(uint8_t c) {
    return RE_map[(uint8_t)c];
}

static const char *RE_tags[] = {
    "U", "S", "E", "C",
    "R",
    "N",
    "X",
    "A", "B", "V", "I"
};

const char *RE_tag_name(int i) {
    return(RE_tags[i]);
}

static const char *exon_type_names[] = {
    "Unknown",
    "Splice", "Exon", "ExonIntron", 
    "IntronRetain",
    "Intron",
    "Exclude",
    "Antisense",
    "AntisenseIntron",
    "Ambiguous",
    "Intergenic"
};

const char *exon_type_name(int i) {
    return(exon_type_names[i]);
}
