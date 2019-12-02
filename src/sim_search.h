// Simimarity Search of short DNA sequences 10~16bp
#ifndef SIM_SEARCH_HEADER
#define SIM_SEARCH_HEADER

typedef uint64_t  base64_t;
typedef uint32_t  base32_t;

typedef struct similarity_search_aux ss_t;

extern ss_t *ss_init();
extern char *ss_query(ss_t *S, char *seq, int e, int *i);
extern int ss_push(ss_t *S, char *seq);
extern void ss_destroy(ss_t *);

#endif
    


   
