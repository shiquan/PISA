/*  anno_thread_pool.c -- A pool of generic worker threads

The original thread_pool.c was copyrighted :

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
// Author : shiquan@genomics.cn

#include <signal.h>
#include <pthread.h>
#include <sys/time.h>
#include <unistd.h>
#include <limits.h>
#include "thread_pool.h"

static int thread_pool_add_result(struct thread_pool_job *j, void *data);
static void wake_next_worker(struct thread_pool_process *q, int locked);
static struct thread_pool_result *thread_pool_next_result_locked(struct thread_pool_process *q);
static void *thread_pool_worker(void *arg);


// A process-queue to hold results from the thread pool.
//
// Each thread pool may have jobs of multiple types being queued up and interleaved, so
// we attach several job process-queues to a single pool.
//
// The jobs themseleves are expected to push their results onto their appropriate
// results queue.

// Add a result to the end of the process result queue.
//
// Returns 0 on success;
//        -1 on failure
static int thread_pool_add_result(struct thread_pool_job *j, void *data) {
    struct thread_pool_process *q = j->q;
    struct thread_pool_result *r;

    pthread_mutex_lock(&q->p->pool_mutex);

    if ( --q->n_processing == 0 )
        pthread_cond_signal(&q->non_processing_c);

    // No results queue is fine if we don't want any results back
    if ( q->in_only ) {
        pthread_mutex_unlock(&q->p->pool_mutex);
        return 0;
    }

    if ( !(r = malloc(sizeof(*r))) )
        return -1;

    r->next = NULL;
    r->data = data;
    r->serial = j->serial;

    q->n_output++;
    if ( q->output_tail ) {
        q->output_tail->next = r;
        q->output_tail = r;
    }
    else {
        q->output_head = q->output_tail = r;
    }

    pthread_cond_broadcast(&q->output_avail_c);

    pthread_mutex_unlock(&q->p->pool_mutex);

    return 0;
}

static struct thread_pool_result *thread_pool_next_result_locked(struct thread_pool_process *q) {
    struct thread_pool_result *r, *last;

    if ( q->shutdown )
        return NULL;

    for ( last = NULL, r = q->output_head; r; last = r, r = r->next) {
        if ( r->serial == q->next_serial )
            break;
    }

    if ( r ) {
        // Remove r from out linked list.
        if ( q->output_head == r )
            q->output_head = r->next;
        else
            last->next = r->next;

        if ( q->output_tail == r )
            q->output_tail = last;

        if ( !q->output_head )
            q->output_tail = NULL;

        q->next_serial++;
        q->n_output--;

        if ( q->qsize && q->n_output < q->qsize ) {
            // Not technically input fill, but can guarantee there is room form the
            // input to go somewhere so we still signal. The waiting code will then
            // check the condition again.
            pthread_cond_signal(&q->input_not_full_c);
            if ( !q->shutdown )
                wake_next_worker(q, 1);
        }
    }
    return r;
}

// Pulls the next item off the process resuly queue. The caller should free it (and
// any internals as appropriate) after use. This doesn't wait for a result to be
// present.
//
// Results will be returned in strict order.
//
// Results thread_pool_result pointer if a result is ready.
//         NULL if not.
struct thread_pool_result *thread_pool_next_result(struct thread_pool_process *q) {
    struct thread_pool_result *r;

    pthread_mutex_lock(&q->p->pool_mutex);
    r = thread_pool_next_result_locked(q);
    pthread_mutex_unlock(&q->p->pool_mutex);
    return r;
}

// Pulls the next item off the process result queue. The caller should free it (and
// any internals as appropriate) after use. This will wait for a result to be present 
// if none are currently available.
//
// Results will be returned in strict order.
//
// Returns thread_pool_resuly pointer if a result is ready.
//         NULL on error or during shutdown.
struct thread_pool_result *thread_pool_next_result_wait(struct thread_pool_process *q) {
    struct thread_pool_result *r;

    pthread_mutex_lock(&q->p->pool_mutex);

    while ( !(r = thread_pool_next_result_locked(q)) ) {
        struct timeval now;
        struct timespec timeout;
        gettimeofday(&now, NULL);
        timeout.tv_sec = now.tv_sec + 10;
        timeout.tv_nsec = now.tv_usec * 1000;

        q->ref_count++;
        if ( q->shutdown ) {
            int rc = --q->ref_count;
            pthread_mutex_unlock(&q->p->pool_mutex);
            if ( rc == 0 )
                thread_pool_process_destroy(q);
            return NULL;
        }
        pthread_cond_timedwait(&q->output_avail_c, &q->p->pool_mutex, &timeout);

        q->ref_count--;
    }
    pthread_mutex_unlock(&q->p->pool_mutex);

    return r;
}

// Returns true if there are no items in the process results queue and also none still
// pending.
int thread_pool_process_empty(struct thread_pool_process *q) {
    int empty;
    pthread_mutex_lock(&q->p->pool_mutex);
    empty = q->n_input == 0 && q->n_processing == 0 && q->n_output == 0;
    pthread_mutex_unlock(&q->p->pool_mutex);
    return empty;
}

void thread_pool_process_ref_incr(struct thread_pool_process *q) {
    pthread_mutex_lock(&q->p->pool_mutex);
    q->ref_count++;
    pthread_mutex_unlock(&q->p->pool_mutex);
}

void thread_pool_process_ref_decr(struct thread_pool_process *q) {
    pthread_mutex_lock(&q->p->pool_mutex);
    if ( --q->ref_count <= 0 ) {
        pthread_mutex_unlock(&q->p->pool_mutex);
        thread_pool_process_destroy(q);
        return;
    }

    // maybe also call destroy here if needed?
    pthread_mutex_unlock(&q->p->pool_mutex);
}

// Returns the number of completed jobs in the process results queue.
int thread_pool_process_len(struct thread_pool_process *q) {
    int len;
    pthread_mutex_lock(&q->p->pool_mutex);
    len = q->n_output;
    pthread_mutex_unlock(&q->p->pool_mutex);

    return len;
}

// Returns the number of completed jobs in the process results queue plus the number
// running and queued up to run.
int thread_pool_process_sz(struct thread_pool_process *q) {
    int len;

    pthread_mutex_lock(&q->p->pool_mutex);
    len = q->n_output + q->n_input + q->n_processing;
    pthread_mutex_unlock(&q->p->pool_mutex);
    return len;
}

// Shutdown a process.
//
// This sets the shutdown flag and wakes any threads waiting on process condition
// variables.
void thread_pool_process_shutdown(struct thread_pool_process *q) {
    pthread_mutex_lock(&q->p->pool_mutex);
    q->shutdown = 1;
    pthread_cond_broadcast(&q->output_avail_c);
    pthread_cond_broadcast(&q->input_not_full_c);
    pthread_cond_broadcast(&q->input_empty_c);
    pthread_cond_broadcast(&q->non_processing_c);
    pthread_mutex_unlock(&q->p->pool_mutex);
}

// Free a result 'r' and if free_data is true also frees the internal r->data result too
void thread_pool_delete_result(struct thread_pool_result *r, int free_data) {
    if ( !r )
        return;

    if ( free_data && r->data)
        free(r->data);

    free(r);
}

// Return the data portion of a thread_pool_result, corresponding to the actual
// "result" itself.
void *thread_pool_result_data(struct thread_pool_result *r) {
    return r->data;
}

// Initialises a thread process-queue.
//
// In_only, if true, indicated that the process generates does not need to hold any
// output. Otherwise an output queue is used to store the results of processing each
// input job.
//
// Returns thread_pool_process pointer on success;
//         NULL on failure.
struct thread_pool_process *thread_pool_process_init(struct thread_pool *p, int qsize,
                                                     int in_only) {
    struct thread_pool_process *q = malloc(sizeof(*q));

    pthread_cond_init(&q->output_avail_c, NULL);
    pthread_cond_init(&q->input_not_full_c, NULL);
    pthread_cond_init(&q->input_empty_c, NULL);
    pthread_cond_init(&q->non_processing_c, NULL);

    q->p = p;
    q->input_head = NULL;
    q->input_tail = NULL;
    q->output_head = NULL;
    q->output_tail = NULL;
    q->next_serial = 0;
    q->curr_serial = 0;
    q->n_input = 0;
    q->n_output = 0;
    q->n_processing = 0;
    q->qsize = qsize;
    q->in_only = in_only;
    q->shutdown = 0;
    q->wake_dispatch = 0;
    q->ref_count = 1;
    q->next = NULL;
    q->prev = NULL;

    thread_pool_process_attach(p, q);
    return q;
}

// Deallocates memeory for a thread process-queue.
// Must be called before the thread pool is destroyed.
void thread_pool_process_destroy(struct thread_pool_process *q) {
    if ( !q )
        return;

    // Ensure it's fully drained before destroying the queue.
    thread_pool_process_reset(q, 0);
    pthread_mutex_lock(&q->p->pool_mutex);
    thread_pool_process_detach(q->p, q);
    thread_pool_process_shutdown(q);

    // Maybe a worker is scanning this queue, so delay destruction.
    if ( --q->ref_count > 0 ) {
        pthread_mutex_unlock(&q->p->pool_mutex);
        return;
    }

    pthread_cond_destroy(&q->output_avail_c);
    pthread_cond_destroy(&q->input_not_full_c);
    pthread_cond_destroy(&q->input_empty_c);
    pthread_cond_destroy(&q->non_processing_c);
    pthread_mutex_unlock(&q->p->pool_mutex);

    free(q);
}

// Attach and detach a thread process-queue with / from the thread pool scheduler.
//
// We need to do attach after making a thread process, but may also wish to temporarily
// detach if we wish to stop running jobs on a specific process while permitting other
// process to continue.
void thread_pool_process_attach(struct thread_pool *p, struct thread_pool_process *q) {
    pthread_mutex_lock(&p->pool_mutex);
    if ( p->q_head ) {
        q->next = p->q_head;
        q->prev = p->q_head->prev;
        p->q_head->prev->next = q;
        p->q_head->prev = q;
    }
    else {
        q->next = q;
        q->prev = q;
    }
    p->q_head = q;
    assert(p->q_head && p->q_head->prev && p->q_head->next);
    pthread_mutex_unlock(&p->pool_mutex);
}
void thread_pool_process_detach(struct thread_pool *p, struct thread_pool_process *q) {
    pthread_mutex_lock(&p->pool_mutex);
    if ( !p->q_head || !q->prev || !q->next)
        goto done;

    struct thread_pool_process *curr = p->q_head, *first = curr;
    do {
        if ( curr == q ) {
            q->next->prev = q->prev;
            q->prev->next = q->next;
            p->q_head = q->prev = NULL;

            // Last one
            if ( p->q_head == q )
                p->q_head = NULL;
            break;
        }
        curr = curr->next;
        
    } while (curr != first);

  done:
    pthread_mutex_unlock(&p->pool_mutex);
}

// A worker thread.
//
// Once woken, each thread checks each process-queue in the pool in turn, looking for
// input jobs that also have room for the output (if it requires storing). If found,
// we execute it and repeat.
// If we checked all input queues and find no such jobs, then we wait until we are
// signalled to check again.
//
static void *thread_pool_worker(void *arg) {

    struct thread_pool_worker *w = (struct thread_pool_worker *)arg;
    struct thread_pool *p = w->p;
    struct thread_pool_job *j;

    for ( ;; ) {
        // Pop an item off the pool queue
        pthread_mutex_lock(&p->pool_mutex);

        assert(p->q_head == 0 || (p->q_head->prev && p->q_head->next));

        int work_to_do = 0;
        struct thread_pool_process *first = p->q_head, *q = first;

        do {
            if ( p->shutdown )
                break;

            // Iterate over queues,
            // finding one with jobs and also room to put the result.          
            if ( q && q->input_head
                 && q->qsize - q->n_output > p->tsize - p->n_waiting) {
                work_to_do = 1;
                break;
            }

            if ( q ) q = q->next;
        } while (q && q != first);

        if ( p->shutdown ) {
          shutdown:
            pthread_mutex_unlock(&p->pool_mutex);
            pthread_exit(NULL);
        }

        if ( !work_to_do ) {
            // We scanned all queues and cannot process any, so we wait.
            p->n_waiting++;

            // Push this thread to the top of the waiting stack.
            if ( p->t_stack_top == -1 || p->t_stack_top > w->idx)
                p->t_stack_top = w->idx;

            p->t_stack[w->idx] = 1;

            pthread_cond_wait(&w->pending_c, &p->pool_mutex);
            p->t_stack[w->idx] = 0;

            // Find new t_stack_top
            int i;
            p->t_stack_top = -1;
            for ( i = 0; i < p->tsize; i++ ) {
                if ( p->t_stack[i]) {
                    p->t_stack_top = i;
                    break;
                }
            }

            p->n_waiting--;
            pthread_mutex_unlock(&p->pool_mutex);
            continue;
        }

        // Otherwise work to do, so process as many items in this queue as possible
        // before switching to another queue. This means threads often end up being
        // dedicated to one type of work.
        q->ref_count++;
        while ( q->input_head && q->qsize - q->n_output > q->n_processing ) {
            if ( p->shutdown )
                goto shutdown;

            j = q->input_head;
            assert(j->p == p);

            if ( !(q->input_head = j->next) )
                q->input_tail = NULL;

            // Transitioning from fill queue to not-full means we can wake up any
            // blocked dispatch threads. We broadcast this as it's only happening once
            // (on the transition) rather than every time we are below qsize.
            q->n_processing++;
            if ( q->n_input-- >= q->qsize )
                pthread_cond_broadcast(&q->input_not_full_c);

            if ( q->n_input == 0 )
                pthread_cond_signal(&q->input_empty_c);

            // Total number of jobs; used to adjust to CPU scaling.
            p->n_jobs--;

            pthread_mutex_unlock(&p->pool_mutex);

            //
            thread_pool_add_result(j, j->func(j->arg, w->idx));
            free(j);

            pthread_mutex_lock(&p->pool_mutex);            
        }

        if ( --q->ref_count == 0 )
            thread_pool_process_destroy(q);
        else
            // Out of jobs on this queue, so restart search from next one.
            // This is equivalent to "work stead"
            p->q_head = q->next;

        pthread_mutex_unlock(&p->pool_mutex);
    }
    return NULL;
}
static void wake_next_worker(struct thread_pool_process *q, int locked) {
    struct thread_pool *p = q->p;
    if ( !locked )
        pthread_mutex_lock(&p->pool_mutex);

    // Update the q_head to be this queue so we'll start processing the queue
    // we know to have results.
    // attached
    assert(q->prev && q->next);
    p->q_head = q;

    // Wake up if we have more jobs waiting than CPUs. This partially combats CPU
    // frequency scaling effects. Starting too many threads and then running out
    // of jobs can cause each thread to have lots of start/stop cycles, which then
    // translates often to CPU frequency scaling adjustments. Instead it is better
    // to only start as many threads as we need to keep the throughput up, meaning
    // some threads run flat out and others are idle.
    //
    // This isn't perfect as we need to know how many can actually start, rather
    // than how many are waiting. A limit on output queue size makes these two
    // figures different.
    assert(p->n_jobs >= q->n_input);

    //int running = p->tsize - p->n_waiting;
    int sig = p->t_stack_top >= 0 && p->n_jobs > p->tsize - p->n_waiting
        && (!q || q->n_processing < q->qsize - q->n_output);

#ifdef AVG_USAGE
    // Track average number of running threads and try to keep close. We permit this
    // to change, but slowly. This avoids "boom and bust" cycles where we read a lot
    // of data, start a lot of jobs, then become idle again. This way some threads run
    // steadily and others dormant, which is better for throughput.
    //
    // It's 50:50 if this is a good thing. It helps some tasks quite significantly
    // while slightlt hindering other (perhaps more usual) jobs.
    if ( ++p->n_count == 256 ) {
        p->n_count >>= 1;
        p->n_running >>= 1;        
    }
    p->n_running += running;
    // Built in lag to avoid see-sawing. Is this safe in all cases?
    if ( sig && p->n_count >= 128 && running * p->n_count > p->n_running+1) sig = 0;
#endif

    if ( sig ) 
        pthread_cond_signal(&p->t[p->t_stack_top].pending_c);

    if ( !locked )
        pthread_mutex_unlock(&p->pool_mutex);
}

// Create a worker pool with n worker threads.
// Returns pool pointer on success;
// NUll on failure.
struct thread_pool *thread_pool_init(int n_threads)
{
    int i;
    struct thread_pool *p = malloc(sizeof(*p));
    p->tsize = n_threads;
    p->n_jobs = 0;
    p->n_waiting = 0;
    p->shutdown = 0;
    p->q_head = NULL;
    p->t_stack = NULL;
    p->n_count = 0;
    p->n_running = 0;
    p->t = malloc(n_threads*sizeof(p->t[0]));

    pthread_mutexattr_t attr;
    pthread_mutexattr_init(&attr);
    pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_RECURSIVE);
    pthread_mutex_init(&p->pool_mutex, &attr);
    pthread_mutexattr_destroy(&attr);

    if ( !(p->t_stack = malloc(n_threads*sizeof(*p->t_stack))) )
        return NULL;

    p->t_stack_top = -1;
    pthread_mutex_lock(&p->pool_mutex);

    for ( i = 0; i < n_threads; i++ ) {
        struct thread_pool_worker *w = &p->t[i];
        p->t_stack[i] = 0;
        w->p = p;
        w->idx = i;
        pthread_cond_init(&w->pending_c, NULL);
        if ( 0 != pthread_create(&w->tid, NULL, thread_pool_worker, w)) {
            pthread_mutex_unlock(&p->pool_mutex);
            return NULL;
        }
    }
    pthread_mutex_unlock(&p->pool_mutex);    
    return p;
}

// Returns the number of requested threads for a pool.
int thread_pool_size(struct thread_pool *p) {
    return p->tsize;
}

// Adds an item to the work pool.
//
// Returns 0 on success
//        -1 on failure
int thread_pool_dispatch(struct thread_pool *p, struct thread_pool_process *q,
                         void *(*func)(void *arg, int idx), void *arg) {
    return thread_pool_dispatch2(p, q, func, arg, 0);
}

// As above but optional non-block flag.
//
// nonblock  0 => block if input queue if full
// nonblock +1 => don't block if input queue if full, but do not add task
// nonblock -1 => add task regardless of whether queue if full (over-size)
int thread_pool_dispatch2(struct thread_pool *p, struct thread_pool_process *q,
                          void *(*func)(void *arg, int idx), void *arg, int nonblock) {
    
    struct thread_pool_job *j;
    pthread_mutex_lock(&p->pool_mutex);
    
    if ( q->n_input >= q->qsize && nonblock == 1 ) {
        pthread_mutex_unlock(&p->pool_mutex);
        errno = EAGAIN;
        return -1;
    }
    
    if ( !(j = malloc(sizeof(*j))) ) {
        pthread_mutex_unlock(&p->pool_mutex);
        return -1;
    }
    j->func = func;
    j->arg = arg;
    j->next = NULL;
    j->p = p;
    j->q = q;
    j->serial = q->curr_serial++;

    if ( nonblock == 0 ) {
        while ( q->n_input >= q->qsize && !q->shutdown && !q->wake_dispatch )
            pthread_cond_wait(&q->input_not_full_c, &q->p->pool_mutex);

        if ( q->shutdown ) {
            free(j);
            pthread_mutex_unlock(&p->pool_mutex);
            return -1;
        }

        if ( q->wake_dispatch ) {
            q->wake_dispatch = 0;
        }
    }

    // total across all queues
    p->n_jobs ++;
    // queue specific
    q->n_input++;

    if ( q->input_tail ) {
        q->input_tail->next = j;
        q->input_tail = j;        
    }
    else {
        q->input_head = q->input_tail = j;
    }

    // Let a worker know we have data.
    //
    // Keep incoming queue at 1 per running thread, so there is always something waiting
    // when they end their current task. If we go above this signal to start more
    // threads (if available). This has the effect of concentrating jobs to fewer cores
    // when we are I/O bound, which in turn benefits systems with auto CPU frequency
    // scaling.
    if ( !q->shutdown )
        wake_next_worker(q, 1);

    pthread_mutex_unlock(&p->pool_mutex);

    return 0;
}

// Wakes up a single thread stuck in dispatch and make it return with error EAGAIN
void thread_pool_wake_dispatch(struct thread_pool_process *q) {
    pthread_mutex_lock(&q->p->pool_mutex);
    q->wake_dispatch = 1;
    pthread_cond_signal(&q->input_not_full_c);
    pthread_mutex_unlock(&q->p->pool_mutex);
}

// Flushes the process-queue, but doesn't exit. This simply drains the queue and ensures
// all worker threads have finished their current tasks associated with this process.
//
// NOT: This does not mean the worker threads are not executing jobs in another process-
// queue.
//
// Return 0 on success;
//       -1 on failure.
int thread_pool_process_flush(struct thread_pool_process *q) {
    int i;
    struct thread_pool *p = q->p;

    // Drains the queue
    pthread_mutex_lock(&p->pool_mutex);

    // Wake up everything for the final sprint!
    for ( i = 0; i < p->tsize; i++ ) {
        if ( p->t_stack[i] )
            pthread_cond_signal(&p->t[i].pending_c);
    }

    // Ensure these is room for the final sprint.
    // Shouldn't be possible to get here, but just incase.
    if ( q->qsize < q->n_output + q->n_input + q->n_processing)
        q->qsize = q->n_output + q->n_input + q->n_processing;

    // Wait for n_input and n_processing to hit zero.
    while (q->n_input || q->n_processing) {
        while (q->n_input )
            pthread_cond_wait(&q->input_empty_c, &p->pool_mutex);
        if ( q->shutdown )
            break;
        while (q->n_processing )
            pthread_cond_wait(&q->non_processing_c, &p->pool_mutex);
        if ( q->shutdown )
            break;
    }

    pthread_mutex_unlock(&p->pool_mutex);
    return 0;
}

// Reset a process to the intial state.
//
// This removes any queued up input jobs, disables any notification of new results/
// output, flushes what is left and then discards any queued output. Anything consumer
// stuck in a wait on results to appear should stay stuck and will only wake up when
// new data is pushed through the queue.
//
// Returns 0 on success;
//        -1 on failure
int thread_pool_process_reset(struct thread_pool_process *q, int free_results) {
    pthread_mutex_lock(&q->p->pool_mutex);
    // prevent next_result from returning data during our flush
    q->next_serial = INT_MAX;

    // Purge any queued input not yet being acted upon
    struct thread_pool_job *j, *jn;
    for ( j = q->input_head; j; j = jn) {
        jn = j->next;
        free(j);
    }
    q->input_head = q->input_tail = NULL;
    q->n_input = 0;

    // Purge any queued output, thus ensuring we have room to flush
    struct thread_pool_result *r, *rn;
    for ( r = q->output_head; r; r = rn) {
        rn = r->next;
        thread_pool_delete_result(r, free_results);
    }
    q->output_head = q->output_tail = NULL;
    q->n_output = 0;
    pthread_mutex_unlock(&q->p->pool_mutex);

    // Wait for any jobs being processed to complete.
    // (TODO: consider how to cancel any currently processing jobs. Probably this is
    // too hard.)
    if ( thread_pool_process_flush(q) != 0 )
        return -1;

    // Discard any new output.
    pthread_mutex_lock(&q->p->pool_mutex);
    for ( r = q->output_head; r; r = rn) {
        rn = r->next;
        thread_pool_delete_result(r, free_results);
    }
    q->output_head = q->output_tail = NULL;
    q->n_output = 0;

    // Finally reset the serial back to the starting point.
    q->next_serial = q->curr_serial = 0;
    pthread_cond_signal(&q->input_not_full_c);
    pthread_mutex_unlock(&q->p->pool_mutex);

    return 0;
}

int thread_pool_process_qsize(struct thread_pool_process *q) {
    return q->qsize;
}

// Destroys a thread pool. The threads are joined into the main thread so they will
// finish their current work load.
void thread_pool_destroy(struct thread_pool *p) {
    int i;
    // Send shutdown message to worker threads
    pthread_mutex_lock(&p->pool_mutex);
    p->shutdown = 1;

    for ( i = 0; i < p->tsize; i++ )
        pthread_cond_signal(&p->t[i].pending_c);

    pthread_mutex_unlock(&p->pool_mutex);

    for ( i = 0; i < p->tsize; i++ )
        pthread_join(p->t[i].tid, NULL);

    pthread_mutex_destroy(&p->pool_mutex);
    for ( i = 0; i < p->tsize; i++ )
        pthread_cond_destroy(&p->t[i].pending_c);

    if ( p->t_stack )
        free(p->t_stack);

    free(p->t);
    free(p);
}

// Destroy a thread pool without waiting on jobs to complete. Use thread_pool_kill(p) to
// quickly exit after a fatal error.
void thread_pool_kill(struct thread_pool *p) {
    int i;
    for ( i = 0; i < p->tsize; i++ ) 
        pthread_kill(p->t[i].tid, SIGINT);

    pthread_mutex_destroy(&p->pool_mutex);
    for ( i = 0; i < p->tsize; i++ )
        pthread_cond_destroy(&p->t[i].pending_c);

    if ( p->t_stack )
        free(p->t_stack);

    free(p->t);
    free(p);    
}

#ifdef THREAD_MAIN_TEST

void *doit(void *arg, int idx)
{
    fprintf(stdout, "idx : %d\n", idx);
    fflush(stdout);
    
    int job = *(int*)arg;
    int *res;
    usleep(500000*((job&31)==31) + random() % 10000);
    res = malloc(sizeof(*res));
    *res = (job<0)? -job : job;
    free(arg);
    return res;
}

void *doit_unorder(void *arg, int idx) {
    fprintf(stdout, "idx : %d\n", idx);
    fflush(stdout);
    int job = *(int*)arg;
    printf("RESULT : %d\n", job);
    return NULL;
}

int test_func_unorder(int n ) {
    struct thread_pool *p = thread_pool_init(n);
    struct thread_pool_process *q = thread_pool_process_init(p, n*2, 1);
    int i;
    for ( i = 0; i < 10000; i++) {
        int *ip = malloc(sizeof(*ip));
        *ip = i;
        thread_pool_dispatch(p, q, doit_unorder, ip);
    }

    thread_pool_process_flush(q);
    thread_pool_process_destroy(q);
    thread_pool_destroy(p);
    return 0;
}

int test_func(int n) {
    int *d = malloc(sizeof(*d));
    *d = 0;    
    struct thread_pool *p = thread_pool_init(n);
    struct thread_pool_process *q = thread_pool_process_init(p, n*2, 0);
    struct thread_pool_result *r;
    int i;
    for ( i = 0; i < 1000; i++ ) {
        int *ip = malloc(sizeof(*ip));
        *ip = i;
        int blk;
        do {
            blk = thread_pool_dispatch2(p, q, doit, ip, 1);
            if ((r = thread_pool_next_result(q))) {
                printf("%d\n", *(int*)r->data);
                thread_pool_delete_result(r, 1);
            }
            if ( blk == -1 ) {
                putchar('.'); fflush(stdout);
                usleep(10000);
            }
        } while (blk == -1);
    }

    thread_pool_process_flush(q);

    while (( r = thread_pool_next_result(q) )) {
        printf("%d\n", *(int*)r->data);
        thread_pool_delete_result(r, 1);
    }

    thread_pool_process_destroy(q);
    thread_pool_destroy(p);
    return 0;
}

int main(int argc, char **argv)
{
    if ( argc != 2 )
        error("thread_pool n_thread");

    int n = atoi(argv[1]);
    if ( n < 1 ) n = 1;
    return test_func_unorder(n);
}

#endif
