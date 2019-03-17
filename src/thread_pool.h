/*
The original thread_pool.h was copyrighted :

    Copyright (c) 2013-2017 Genome Research Ltd.

    Author: James Bonfield <jkb@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

// This file was edited to access worker id (thread layer) in work data (process layer).
// The struct names have been changed to avoid potential conflict with JB's thread_pool.

#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include "utils.h"
#include <pthread.h>
#include <stdint.h>

struct thread_pool; 
struct thread_pool_process; 
struct thread_pool_job;
struct thread_pool_result; 
struct thread_pool_worker;

struct thread_pool_job {
    void *(*func)(void *arg, int idx);
    void *arg;
    struct thread_pool_job *next;
    struct thread_pool     *p;
    struct thread_pool_process *q;
    uint64_t serial;
};

struct thread_pool_result {
    struct thread_pool_result *next;
    // sequential number for ordering
    uint64_t serial;
    // result itself
    void *data; 
};

struct thread_pool_worker {
    struct thread_pool *p;
    int idx;
    pthread_t tid;
    // when waiting for a job
    pthread_cond_t pending_c;
};

// An IO queue consists of a queue of jobs of execute (the "input" side) and a queue
// of job results post-execution (the "output" side).
// We have size limits to prevent either queue from growing too large and serial numbers
// to ensure sequential consumption of the output.
// The thread pool may have many hetergeneous tasks, each using its own io_queue mixed
// into the sam thread pool.
struct thread_pool_process {
    struct thread_pool *p;
    // input list
    struct thread_pool_job     *input_head;
    struct thread_pool_job     *input_tail;
    // output list
    struct thread_pool_result  *output_head;
    struct thread_pool_result  *output_tail;
    // max size of I/O queues
    int qsize;
    // next serial for output
    uint64_t next_serial;
    // current serial (next input)
    uint64_t curr_serial;

    // no. items in input queue; was njobs
    int n_input;
    // no. items in output queue
    int n_output;
    // no. ithems being processed (excuting)
    int n_processing;

    // true if pool is being destroyed
    int shutdown;
    // if true, don't queue result up.
    int in_only;
    // unblocks waiting dispatchers
    int wake_dispatch;

    // used to track safe destruction
    int ref_count;

    // signalled on each new output
    pthread_cond_t output_avail_c;
    // input queue is no longer full
    pthread_cond_t input_not_full_c;
    // input queue has become empty
    pthread_cond_t input_empty_c;
    // n_processing has hit zero
    pthread_cond_t non_processing_c;

    // to form circular linked list
    struct thread_pool_process *next, *prev;    
};


// The single pool structure itself.
//
// This knows nothing about the nature of the jobs or where their output is going,
// but it maintains a list of queues associated with this pool from which the jobs
// are taken.
struct thread_pool {
    // how many workers waiting for new jobs
    int n_waiting;
    // how many total jobs are waiting in all queues
    int n_jobs;
    // ture if pool is be destroyed
    int shutdown;

    // I/O queues to check for jobs in and to put results.
    // Forms a circular linked list. (q_head may be amended to
    // point to the most recently updated.)
    struct thread_pool_process *q_head;

    // threads
    // maxinum number of jobs
    int tsize;
    struct thread_pool_worker *t;
    // array of worker IDs free
    int *t_stack, t_stack_top;

    // A single mutex used when updateing this and any associated structure.
    pthread_mutex_t pool_mutex;

    // Tracking of average number of running jobs.
    // This can be used to dampen any hysteresis caused by bursty input availability.
    int n_count, n_running;

};

struct thread_pool_result *thread_pool_next_result(struct thread_pool_process *q);
struct thread_pool_result *thread_pool_next_result_wait(struct thread_pool_process *q);
int thread_pool_process_empty(struct thread_pool_process *q);
void thread_pool_process_ref_incr(struct thread_pool_process *q);
void thread_pool_process_ref_decr(struct thread_pool_process *q);
int thread_pool_process_len(struct thread_pool_process *q);
int thread_pool_process_sz(struct thread_pool_process *q);
void thread_pool_process_shutdown(struct thread_pool_process *q);
void thread_pool_delete_result(struct thread_pool_result *, int);
void *thread_pool_result_data(struct thread_pool_result *);
struct thread_pool_process *thread_pool_process_init(struct thread_pool *, int, int);
void thread_pool_process_destroy(struct thread_pool_process *q);
void thread_pool_process_attach(struct thread_pool *p, struct thread_pool_process *q);
void thread_pool_process_detach(struct thread_pool *p, struct thread_pool_process *q);
struct thread_pool *thread_pool_init(int n_threads);
int thread_pool_size(struct thread_pool *p);
int thread_pool_dispatch(struct thread_pool *p, struct thread_pool_process *q,void *(*func)(void *arg, int idx),void *arg);
int thread_pool_dispatch2(struct thread_pool *p, struct thread_pool_process *q,void *(*func)(void *arg, int idx), void *arg, int nonblock);
void thread_pool_wake_dispatch(struct thread_pool_process *q);
int thread_pool_process_flush(struct thread_pool_process *q);
int thread_pool_process_reset(struct thread_pool_process *q, int free_results);
int thread_pool_process_qsize(struct thread_pool_process *q);
void thread_pool_destroy(struct thread_pool *p);
void thread_pool_kill(struct thread_pool *p);

#endif
