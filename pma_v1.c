#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <pthread.h>
#define NO_OF_THREADS 10
/* Reserve 8 bits to allow for fixed point arithmetic. */
#define MAX_SIZE ((1ULL << 56) - 1ULL)
/* Height-based (as opposed to depth-based) thresholds. */
/* Upper density thresholds. */
static const double t_h = 0.75; /* root. */
static const double t_0 = 1.00; /* leaves. */
/* Lower density thresholds. */
static const double p_h = 0.20; /* root. */
static const double p_0 = 0.25; /* leaves. */

static const uint8_t max_sparseness = 1 / p_0;
static const uint8_t largest_empty_segment = 1 * max_sparseness;

struct _pma;
typedef struct _pma pma_t, *PMA;
/*
 * Returns the 1-based index of the last (i.e., most significant) bit set in x.
 */
uint64_t last_bit_set(uint64_t x)
{
    // assert (x > 0);
    return (sizeof(uint64_t) * 8 - __builtin_clzll(x));
}

uint64_t floor_lg(uint64_t x)
{
    return (last_bit_set(x) - 1);
}

uint64_t ceil_lg(uint64_t x)
{
    return (last_bit_set(x - 1));
}

/* Returns the largest power of 2 not greater than x
 * (i.e., $2^{\lfloor \lg x \rfloor}$).
 */
uint64_t hyperfloor(uint64_t x)
{
    return (1 << floor_lg(x));
}

/* Returns the smallest power of 2 not less than x
 * (i.e., $2^{\lceil \lg x \rceil}$).
 */
uint64_t hyperceil(uint64_t x)
{
    return (1 << ceil_lg(x));
}

uint64_t ceil_div(uint64_t x, uint64_t y)
{
    assert(x > 0);
    return (1 + ((x - 1) / y));
}

typedef uint64_t key_t;
typedef uint64_t val_t;

// struct definitions

typedef struct
{
    int operation;
    uint64_t version;
    key_t key;
    val_t val;
} marker_t;

typedef struct
{
    key_t key;
    val_t val;
    uint64_t version;
    marker_t mark;
} keyval_t;

/* Returns true if keyval is empty and false otherwise. */
bool keyval_empty(const keyval_t *keyval)
{
    return (keyval->key == 0ULL);
}

/* Sets keyval to be empty. */
void keyval_clear(keyval_t *keyval)
{
    keyval->key = 0ULL;
    keyval->val = 0ULL;
}

/* compare and swap implementation for marker*/
bool CASM(marker_t *p, marker_t old, marker_t new)
{
    if ((*p).key == old.key && (*p).val == old.val && (*p).operation == old.operation && (*p).version == old.version)
    {
        *p = new;
        return true;
    }
    else
    {
        return false;
    }
}

/* double compare and swap implementation for value*/
bool DCAS(uint64_t *p, uint64_t *q, uint64_t old_p, uint64_t new_p, uint64_t old_q, uint64_t new_q)
{
    if (*p == old_p && *q == old_q)
    {
        (*p = new_p) && (*q = new_p);
        return true;
    }
    else
    {
        return false;
    }
}

/* compare and swap implementation for value*/
bool CAS(uint64_t *p, uint64_t old, uint64_t new)
{
    if (*p == old)
    {
        *p = new;
        return true;
    }
    else
    {
        return false;
    }
}

PMA pma_create(void);
PMA pma_from_array(keyval_t *array, uint64_t n);
void pma_destroy(PMA *pma);
void insert(PMA pma, key_t key, val_t val, uint64_t idx);
void move(PMA pma, uint64_t from, uint64_t to);
bool rebalance_move(PMA pma, uint64_t from, uint64_t to);
void pma_insert(PMA pma, key_t key, val_t val, uint64_t idx);
void pma_insert_after(PMA pma, int64_t i, key_t key, val_t val);
bool pma_delete(PMA pma, key_t key);
void pma_delete_at(PMA pma, int64_t i);
void pma_get(PMA pma, int64_t i, keyval_t *keyval);
uint64_t pma_capacity(PMA p);
uint64_t pma_count(PMA p);

/* TODO: For testing purposes only. */
uint8_t pma_segment_size(PMA pma);

struct _pma
{
    uint64_t n; /* Number of elements. */
    uint64_t m; /* Size of the array (capacity). */
    uint8_t s;  /* Size of the segments. */
    uint64_t num_segments;
    uint8_t h;      /* Height of the tree. */
    double delta_t; /* Delta for upper density threshold. */
    double delta_p; /* Delta for lower density threshold. */
    keyval_t *array;
};

/* TODO: For testing purposes only. */
uint8_t pma_segment_size(PMA pma)
{
    return (pma->s);
}

static void rebalance(PMA pma, int64_t i);
static bool pack(PMA pma, uint64_t from, uint64_t to, uint64_t n);
static bool spread(PMA pma, uint64_t from, uint64_t to, uint64_t n);
static bool resize(PMA pma);
static uint64_t compute_capacity(PMA pma);

PMA pma_test;

PMA pma_create()
{
    PMA pma = (PMA)malloc(sizeof(pma_t));
    pma->n = 0;
    pma->s = largest_empty_segment;
    pma->m = (10000000);
    pma->num_segments = pma->m / pma->s;
    pma->h = floor_lg(pma->num_segments) + 1;
    pma->delta_t = (t_0 - t_h) / pma->h;
    pma->delta_p = (p_h - p_0) / pma->h;
    pma->array = (keyval_t *)malloc(sizeof(keyval_t) * pma->m);
    for (int x = 0; x < pma->m; x++)
    {
        pma->array[x].key = 0ULL;
        pma->array[x].val = 0ULL;
        marker_t mk = {.operation = 0, .key = 0, .val = 0, .version = 0};
        pma->array[x].mark = mk;
        pma->array[x].version = 0ULL;
    }
    return (pma);
}

void pma_destroy(PMA *pma)
{
    free((*pma)->array);
    (*pma)->array = NULL;
    free(*pma);
    *pma = NULL;
}
void insert(PMA pma, key_t key, val_t val, uint64_t idx)
{
    while (true)
    {
        bool help = false;
        if (pma->array[idx - 1].version < pma->array[idx - 1].mark.version)
        {
            // printf("Version: %d and Marker Version %d\n", pma->array[idx-1].version ,pma->array[idx-1].mark.version);
            continue;
        }
        if (pma->array[idx - 1].mark.operation == 0)
        {
            // perform CAS on the cell marker
            marker_t old = pma->array[idx - 1].mark;
            marker_t new = {.operation = 1, .key = pma->array[idx - 1].key, .val = pma->array[idx - 1].val, .version = old.version + 1};
            if (CASM(&pma->array[idx - 1].mark, old, new))
            {
                // printf("Index: %d\n",idx);
                while (true)
                {
                    uint64_t old_key = pma->array[idx].key;
                    uint64_t old_val = pma->array[idx].val;
                    if (DCAS(&pma->array[idx].key, &pma->array[idx].val, old_key, key, old_val, val))
                    {
                        break;
                    }
                }
                pma->array[idx - 1].mark.operation = 0;
                pma->array[idx - 1].version = pma->array[idx - 1].mark.version;
                break;
            }
            else
            {
                continue;
            }
        }
        // printf("Version: %d and Marker Version %d\n", pma->array[idx-1].version ,pma->array[idx-1].mark.version);
    }
}

uint64_t read(PMA pma, uint64_t idx)
{
    return 0ULL;
}

void move(PMA pma, uint64_t from, uint64_t to)
{
    while (true)
    {
        bool help = false;
        if (pma->array[from].version < pma->array[from].mark.version)
        {
            // printf("Version: %d and Marker Version %d\n", pma->array[i-1].version ,pma->array[i-1].mark.version);
            continue;
        }
        if (pma->array[from].mark.operation == 0)
        {
            // perform CAS on the cell marker
            marker_t old = pma->array[from].mark;
            marker_t new = {.operation = 1, .key = 0ULL, .val = 0ULL, .version = old.version + 1};
            if (CASM(&pma->array[from].mark, old, new))
            {
                while (1)
                {
                    uint64_t old_key = pma->array[to].key;
                    uint64_t old_val = pma->array[to].val;
                    if (DCAS(&pma->array[to].key, &pma->array[to].val, old_key, pma->array[from].key, old_val, pma->array[from].val))
                    {
                        pma->array[to].mark = pma->array[from].mark;
                        pma->array[to].mark.operation = 0;
                        pma->array[to].version = pma->array[from].mark.version;
                        pma->array[from].key = 0ULL;
                        break;
                    }
                }
                pma->array[from].mark.operation = 0;
                pma->array[from].version = pma->array[from].mark.version;
                break;
            }
            else
            {
                continue;
            }
        }
    }
}
bool rebalance_move(PMA pma, uint64_t from, uint64_t to)
{

    // bool help = false;
    // if (pma->array[from].version < pma->array[from].mark.version)
    // {
    //     printf("Version: %d and Marker Version %d\n", pma->array[from].version ,pma->array[from].mark.version);
    //     return false;
    // }
    if (pma->array[from].mark.operation == 0 )
    {
        // perform CAS on the cell marker
        marker_t old = pma->array[from].mark;
        marker_t new = {.operation = 1, .key = 0, .val = 0, .version = old.version + 1};
        if (CASM(&pma->array[from].mark, old, new))
        {
            while (true)
            {
                uint64_t old_key = pma->array[to].key;
                uint64_t old_val = pma->array[to].val;
                if (DCAS(&pma->array[to].key, &pma->array[to].val, old_key, pma->array[from].key, old_val, pma->array[from].val))
                {
                    pma->array[to].mark = pma->array[from].mark;
                    pma->array[to].mark.operation = 0;
                    pma->array[to].version = pma->array[from].mark.version;
                    pma->array[from].key = 0ULL;
                    // printf("Read Index: %d, Write Index: %d, Reader Op: %d, Writer Op: %d\n", from, to, pma->array[from].mark.operation, pma->array[to].mark.operation);
                    break;
                }
            }
            pma->array[from].mark.operation = 0;
            pma->array[from].version = pma->array[from].mark.version;
            // assert(pma->array[from].mark.operation == 0 && pma->array[to].mark.operation == 0);
            return true;
        }
        else
        {
            return false;
        }
    }else{
        return false;
    }
}
void pma_insert(PMA pma, key_t key, val_t val, uint64_t idx)
{
    int64_t j = idx;
    // printf("capacity: %d\n",pma->m);
    while (j < pma->m && !keyval_empty(&(pma->array[j]))){
        j++;
        // printf("Cant insert here!");
        // return;
    }
    
    if (j < pma->m)
    {
        // Found a position to the right
        while (j > idx)
        {
            // printf("Index: %d\n",idx);
            move(pma, j - 1, j);
            j--;
        }
        // printf("Index: %d\n",idx);       
        insert(pma, key, val, idx);
    }
    pma->n++;
    rebalance(pma, idx);
    // printf("Elements: %d\n",pma->n);
    // printf("Current capacity: %d\n",pma->m);
}

static void rebalance(PMA pma, int64_t i)
{
    // printf("HEYY1\n");
    int64_t window_start, window_end;
    uint8_t height = 0;
    uint64_t occupancy = (keyval_empty(&(pma->array[i]))) ? 0 : 1;
    int64_t left_index = i - 1;
    int64_t right_index = i + 1;
    double density, t_height, p_height;
    do
    {
        // printf("Current capacity: %d\n",pma->m);
        uint64_t window_size = pma->s * (1 << height);
        uint64_t window = i / window_size;
        window_start = window * window_size;
        window_end = window_start + window_size;
        while (left_index >= window_start)
        {
            if (!keyval_empty(&(pma->array[left_index])))
                occupancy++;
            left_index--;
        }
        while (right_index < window_end)
        {
            if (!keyval_empty(&(pma->array[right_index])))
                occupancy++;
            right_index++;
        }
        density = (double) occupancy / window_size;
        // printf("Current density, occupancy, window size in this segment: %lf %ld %ld\n",density, occupancy, window_size);
        t_height = t_0 - (height * pma->delta_t);
        p_height = p_0 + (height * pma->delta_p);
        height++;
    } while ((density < p_height ||
              density >= t_height) &&
             height < pma->h);
    /* Found a window within threshold. */
    if (density >= p_height && density < t_height)
    {
        if (!(pack(pma, window_start, window_end, occupancy) && spread(pma, window_start, window_end, occupancy)))
        {
            // printf("PACK & SPREAD FAILED\n");
            rebalance(pma, i);
        }
    }
    else
    {
        printf("Current capacity before resize: %d\n",pma->m);
        if (!resize(pma))
        {
            printf("RESIZE FAILED\n");
            rebalance(pma, i);
        }

    }
}

static bool pack(PMA pma, uint64_t from, uint64_t to, uint64_t n)
{
    assert(from < to);
    uint64_t read_index = from;
    uint64_t write_index = from;
    while (read_index < to)
    {
        if (!keyval_empty(&(pma->array[read_index])))
        {
            // printf("Read Index: %d, Write Index: %d\n", read_index, write_index);
            if (read_index > write_index)
            {
                if (!rebalance_move(pma, read_index, write_index)){
                    // printf("NAH\n");
                    return false;
                }
                // assert(pma->array[write_index].mark.operation == 2 && pma->array[read_index].mark.operation == 0);
            }
            write_index++;
        }
        read_index++;
    }
    return true;
}
static bool spread(PMA pma, uint64_t from, uint64_t to, uint64_t n)
{
    assert(from < to);
    uint64_t capacity = to - from;
    uint64_t frequency = (capacity << 8) / n; /* 8-bit fixed point arithmetic. */
    uint64_t read_index = from + n - 1;
    uint64_t write_index = (to << 8) - frequency;
    while ((write_index >> 8) > read_index)
    {
        // printf("Read Index: %d, Write Index: %d\n", read_index, write_index >> 8);
        if (!rebalance_move(pma, read_index, write_index >> 8))
            return false;
        // assert(pma->array[write_index >> 8].mark.operation == 0 && pma->array[read_index].mark.operation == 0);
        read_index--;
        write_index -= frequency;
    }
    return true;
}

static bool resize(PMA pma)
{
    if (!pack(pma, 0, pma->m, pma->n)){
        printf("PACK FAILED\n");
        return false;
    }
    uint64_t old_m = pma->m;
    uint64_t new_m = compute_capacity(pma);
    // printf("Capacity increased from old: %d, new: %d\n", old_m, new_m);
    pma->h = floor_lg(pma->num_segments) + 1;
    pma->delta_t = (t_0 - t_h) / pma->h;
    pma->delta_p = (p_h - p_0) / pma->h;
    pma->array = (keyval_t *)realloc(pma->array, sizeof(keyval_t) * new_m);
    for (int x = old_m; x < new_m; x++)
    {
        pma->array[x].key = 0;
        pma->array[x].val = 0;
        pma->array[x].version = 0;
        marker_t mk1 = {.operation = 0, .key = 0, .val = 0, .version = 0};
        pma->array[x].mark = mk1;
    }
    pma->m = new_m;
    for (uint64_t x = pma->n; x < pma->m; x++)
    {
        insert(pma, 0, 0, x);
    }
    if (!spread(pma, 0, pma->m, pma->n)){
        printf("SPREAD FAILED\n");
        return false;
    }
    return true;
}

static uint64_t compute_capacity(PMA pma)
{
  pma->s = ceil_lg(pma->n);                     /* Ideal segment size. */
  pma->num_segments = ceil_div(pma->n, pma->s); /* Ideal number of segments. */
  /* The number of segments has to be a power of 2. */
  pma->num_segments = hyperceil(pma->num_segments);
  /* Update the segment size accordingly. */
  pma->s = ceil_div(pma->n, pma->num_segments);
  uint64_t m = pma->s * pma->num_segments;
  /* Scale up as much as possible. */
  m = max_sparseness * pma->m;
  pma->s = max_sparseness * pma->s;
  assert(m <= MAX_SIZE);
  assert(m > pma->n);
  return m;
}

void pma_print(PMA pma)
{
  for (int x = 0; x < pma->m; x++)
  {
    printf("Index: %d, key: %d, value: %d, version: %d, marker_version: %d, marker_operation: %d\n", x, pma->array[x].key, pma->array[x].val, pma->array[x].version, pma->array[x].mark.version, pma->array[x].mark.operation);
  }
}
void *threadFunc(void *args)
{
    int tid = *((int *)args);
    // printf("tid is %d \n", tid);
    for (uint64_t x = (100000)*tid+1; x <= (100000)*(tid+1); x++)
    {
        pma_insert(pma_test, x, x, x);
    }
    pthread_exit(NULL);
}
int main()
{
    pma_test = pma_create();
    pthread_t tid[10];
    int msec = 0, trigger = 10; /* 10ms */
    clock_t before = clock();
    // pma_print(pma_test);
    for (int x = 0; x < NO_OF_THREADS; x++)
    {
        int *temp = calloc(1, sizeof(int));
        *temp = x;
        if (pthread_create(&tid[x], NULL, threadFunc, (void *)temp) != 0)
        {
        printf("Cannot spawn thread\n");
        exit(-1);
        }
    }
    // printf("after pthread_create \n");
    for (int x = 0; x < NO_OF_THREADS; x++)
    {
        if (pthread_join(tid[x], NULL) != 0)
        {
        perror("pthread_join");
        printf("Cannot join thread\n");
        exit(-1);
        }
    }
    clock_t difference = clock() - before;
    msec = difference * 1000 / CLOCKS_PER_SEC;
    printf("Time taken: %d seconds %d milliseconds\n",
           msec / 1000, msec % 1000);
    // pma_print(pma_test);
    printf("Elements: %d\n", pma_test->n);
    printf("Capacity: %d", pma_test->m);
    pthread_exit(NULL);
}