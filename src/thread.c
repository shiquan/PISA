#include "thread.h"
#include "utils.h"


void tpool_thread(void *_p)
{
    struct tpool_work *w;
    struct tpool *p = _p;
    for (;;) {
        pthread_mutex_lock(&p->lock);
        while (p->n == 0 && !p->shutdown)
            pthread_cond_wait(&p->not_empty, &p->lock);
        if (p->shutdown) {
            pthread_mutex_unlock(&p->lock);
            pthread_exit(NULL);
        }
        w = p->head;
        p->n--;
        if (p->n == 0)
            p->head = p->tail = NULL;
        else
            p->head = w->next;

        if (!p->do_not_block && p->n == (p->m-1))
            pthread_cond_broadcast(&p->not_full);
        
        if (p->n == 0)
            pthread_cond_signal(&p->empty);

        pthread_mutex_unlock(&p->lock);
        
        (*(w->routine))(w->arg);
        free(w);
    }
}
struct tpool *tpool_init(int n_thread, int max_size, int do_not_block)
{
    struct tpool *p = malloc(sizeof(*p));
    memset(p, 0, sizeof(*p));
    p->n_thread = n_thread;
    p->m = max_size;
    p->do_not_block = do_not_block;
    p->threads = malloc(sizeof(pthread_t)*n_thread);
    pthread_mutex_init(&p->lock, NULL);
    pthread_cond_init(&p->not_empty, NULL);
    pthread_cond_init(&p->not_full, NULL);
    pthread_cond_init(&p->empty, NULL);
    int i;
    for (i = 0; i < n_thread; ++i)
        pthread_create(&p->threads[i], NULL, tpool_thread, p);

    return p;
}

int tpool_add_work(struct tpool *p, void *routine, void *arg)
{
    struct tpool_work *w;
    pthread_mutex_lock(&p->lock);
    if (p->n == p->m && p->do_not_block) {
        pthread_mutex_unlock(&p->lock);
        return -1;
    }

    while (p->n == p->m && !p->shutdown)
        pthread_cond_wait(&p->not_full, &p->lock);

    if (p->shutdown) {
        pthread_mutex_unlock(&p->lock);
        return -1;
    }

    w = malloc(sizeof(*w));
    w->routine = routine;
    w->arg = arg;
    w->next = NULL;

    if (p->n == 0) {
        p->tail = p->head = w;
        pthread_cond_broadcast(&p->not_empty);
    }
    else {
        p->tail->next = w;
        p->tail = w;
    }
    p->n++;
    pthread_mutex_unlock(&p->lock);
    return 1;
}

void tpool_destroy(struct tpool *p)
{
    pthread_mutex_unlock(&p->lock);

    if (p->shutdown) {
        pthread_mutex_unlock(&p->lock);
        return;
    }

    while (p->n != 0) 
        pthread_cond_wait(&p->empty, &p->lock);
        
    p->shutdown = 1;

    pthread_mutex_unlock(&p->lock);
    pthread_cond_broadcast(&p->not_empty);
    pthread_cond_broadcast(&p->not_full);

    int i;
    for (i = 0; i < p->n_thread; i++)
        pthread_join(p->threads[i], NULL);
    
    free(p->threads);
    free(p);
    return;
}
