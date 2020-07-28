#ifndef REGION_IDX_H
#define REGION_IDX_H


struct region_index;

struct region_itr {
    int n;
    void **rets;
};

struct region_index *region_index_create();
void region_index_destroy(struct region_index *idx);

void index_bin_push(struct region_index *idx, uint32_t start, uint32_t end, void *new);

struct region_itr *region_query(struct region_index *idx, int start, int end);
void region_itr_destroy(struct region_itr *itr);


#endif
