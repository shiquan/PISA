#ifndef PISA_THREAD_H
#define PISA_THREAD_H

#include<stdio.h>
#include<pthread.h>

struct tpool_work {
    void (*routine)();
    void *arg;
    struct tpool_work *next;
};

struct tpool {
    int                n_thread;
    int                m, n; // queue size
    int                do_not_block;
    pthread_t         *threads;
    struct tpool_work *head;
    struct tpool_work *tail;
    pthread_mutex_t    lock;
    pthread_cond_t     not_empty;
    pthread_cond_t     not_full;
    pthread_cond_t     empty;
    int                shutdown;
};

struct tpool *tpool_init(int n_thread, int max_size, int do_not_block);
int tpool_add_work(struct tpool *p, void *routine, void *arg);
void tpool_destroy(struct tpool *p);

#endif
