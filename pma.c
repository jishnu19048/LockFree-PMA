/*
 * Copyright (c) 2014 Pablo Montes <pabmont@gmail.com>
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */
#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <pthread.h>
#define NO_OF_THREADS 2
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

/* Reserve 8 bits to allow for fixed point arithmetic. */
#define MAX_SIZE ((1ULL << 56) - 1ULL)

/* Height-based (as opposed to depth-based) thresholds. */
/* Upper density thresholds. */
static const double t_h = 0.75; /* root. */
static const double t_0 = 1.00; /* leaves. */
/* Lower density thresholds. */
static const double p_h = 0.50; /* root. */
static const double p_0 = 0.25; /* leaves. */

static const uint8_t max_sparseness = 1 / p_0;
static const uint8_t largest_empty_segment = 1 * max_sparseness;

struct _pma;
typedef struct _pma pma_t, *PMA;

PMA pma_create(void);
PMA pma_from_array(keyval_t *array, uint64_t n);
void pma_destroy(PMA *pma);
bool pma_find(PMA pma, key_t key, int64_t *index);
bool pma_insert(PMA pma, key_t key, val_t val);
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

PMA pma_create()
{
  PMA pma = (PMA)malloc(sizeof(pma_t));
  pma->n = 0;
  pma->s = largest_empty_segment;
  pma->m = (1ULL << largest_empty_segment);
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
PMA pma1;

void pma_destroy(PMA *pma)
{
  free((*pma)->array);
  (*pma)->array = NULL;
  free(*pma);
  *pma = NULL;
}

/*
 * Perform a modified binary search, with $O (\lg n)$ comparisons, that allows
 * gaps of size $O(1)$ in the array.
 * Returns true if the element is found and false otherwise.
 * If the element is found, index holds the position in the array where the
 * element associated with key is stored.
 * If the element is not found, index holds the position of the predecessor or
 * -1 if no predecessor exist in the array.
 */
bool pma_find(PMA pma, key_t key, int64_t *index)
{
  int64_t from = 0;
  int64_t to = pma->m - 1;
  while (from <= to)
  {
    int64_t mid = from + (to - from) / 2;
    int64_t i = mid;
    /* Start scanning left until we find a non-empty slot or we reach past the
     * beginning of the subarray. */
    while (i >= from &&
           keyval_empty(&(pma->array[i])))
      i--;
    if (i < from)
    { /* Everything between from and mid (inclusive) is empty. */
      from = mid + 1;
    }
    else
    {
      if (pma->array[i].key == key)
      {
        *index = i;
        return (true);
      }
      else if (pma->array[i].key < key)
      {
        from = mid + 1;
      }
      else
      { /* pma->array [i].key > key */
        to = i - 1;
      }
    }
  }
  /* Didn't find `key'. `to' should hold its predecessor (unless it's empty). */
  *index = to;
  while (*index >= 0 && keyval_empty(&(pma->array[*index])))
    (*index)--;
  return (false);
}

bool pma_insert(PMA pma, key_t key, val_t val)
{
  int64_t i;
  if (!pma_find(pma, key, &i))
  { /* We do not allow duplicates.*/
    // printf("index: %d\n",i);
    if (pma->n == 0)
    {
      // printf("n is 0 and i is %d\n", i + 1);
      while (true)
      {
        if (pma->array[i + 1].mark.operation == 0)
        {
          uint64_t old_key = pma->array[i + 1].key;
          uint64_t old_val = pma->array[i + 1].val;
          if (!DCAS(&pma->array[i + 1].key, &pma->array[i + 1].val, old_key, key, old_val, val))
          {
            continue;
          }
          break;
        }
      }
      pma->n++;
      // printf("1 element added, Version: %d\n", pma->array[i + 1].mark.version);
      return true;
    }
    else
    {
      pma_insert_after(pma, i, key, val);
      return (true);
    }
  }
  return (false);
}

void pma_insert_after(PMA pma, int64_t i, key_t key, val_t val)
{
  assert(i >= -1);
  assert(i < (int64_t)pma->m);
  assert(i >= 0 && !keyval_empty(&(pma->array[i])) || i >= -1);
  /* Find an empty space to the right of i. There should be one close by. */
  int64_t j = i + 1;
  // printf("capacity: %d\n",pma->m);
  while (j < pma->m && !keyval_empty(&(pma->array[j])))
    j++;
  // printf("index: %d\n",i);
  if (j < pma->m)
  { /* Found one. */
    // printf("FOUND\n");
    while (j > i + 1)
    { /* Push elements to make space for the new element. */
      while (true)
      {
        bool help = false;
        if (pma->array[j-1].version < pma->array[j-1].mark.version)
        {
          // printf("Version: %d and Marker Version %d\n", pma->array[j-1].version ,pma->array[j-1].mark.version);
          printf("Stuck1\n");
          continue;
        }
        if (pma->array[j-1].mark.operation == 0)
        {
          marker_t old = pma->array[j-1].mark;
          marker_t new = {.operation = 1, .key = 0, .val = 0, .version = old.version + 1};
          if (CASM(&pma->array[j-1].mark, old, new))
          {
            uint64_t old_key = pma->array[j].key;
            uint64_t old_val = pma->array[j].val;
            if (DCAS(&pma->array[j].key, &pma->array[j].val, old_key, pma->array[j-1].key, old_val, pma->array[j-1].val))
            {
              pma->array[j].mark = pma->array[j-1].mark;
              pma->array[j].mark.operation = 0;
              pma->array[j].version = pma->array[j-1].mark.version;
              uint64_t old_key_from = pma->array[j-1].key;
              help = true;
              if (!CAS(&pma->array[j-1].key, old_key_from, 0ULL))
              {
                // can fail here because can be helped
              }
              (pma->array[j-1].mark.operation = 0);
              (pma->array[j-1].version = pma->array[j-1].mark.version);
              break;
            }
            else
            {
              (pma->array[j-1].mark.operation = 0);
              (pma->array[j-1].version = pma->array[j-1].mark.version);
              continue;
            }
          }
          else
          {
            printf("Stuck2\n");
            continue;
          }
        }
        else
        {
          if (help)
          {
            pma->array[j-1].key = pma->array[j-1].mark.key;
            pma->array[j-1].val = pma->array[j-1].mark.val;
          }
        }
      }
      j--;
    }
    while (true)
    {
      bool help = false;

      if (pma->array[i].version < pma->array[i].mark.version)
      {
        // printf("Version: %d and Marker Version %d\n", pma->array[i].version ,pma->array[i].mark.version);
        printf("Stuck3\n");
        continue;
      }
      if (pma->array[i].mark.operation == 0)
      {
        // perform CAS on the cell marker
        marker_t old = pma->array[i].mark;
        marker_t new = {.operation = 1, .key = pma->array[i].key, .val = pma->array[i].val, .version = old.version + 1};
        if (CASM(&pma->array[i].mark, old, new))
        {
          uint64_t old_key = pma->array[i + 1].key;
          uint64_t old_val = pma->array[i + 1].val;
          if (!DCAS(&pma->array[i + 1].key, &pma->array[i + 1].val, old_key, key, old_val, val))
          {
            pma->array[i].mark.operation = 0;
            pma->array[i].version = pma->array[i].mark.version;
            continue;
          }
          pma->array[i].mark.operation = 0;
          pma->array[i].version = pma->array[i].mark.version;
          break;
        }
        else
        {
          printf("Stuck4\n");
          continue;
        }
      }
      // printf("HELLO\n");
    }
    i = i + 1; /* Update i to point to where the new element is. */
    // printf("verSion: %d and marker version: %d\n", pma->array[i-1].version, pma->array[i-1].mark.version);
  }
  else
  { /* No empty space to the right of i. Try left. */
    j = i - 1;
    while (j >= 0 && !keyval_empty(&(pma->array[j])))
      j--;
    if (j >= 0)
    { /* Found one. */
      while (j < i)
      { /* Push elements to make space for the new element. */
        while (true)
        {
          bool help = false;
          if (pma->array[j+1].version < pma->array[j+1].mark.version)
          {
            // printf("Version: %d and Marker Version %d\n", pma->array[j+1].version ,pma->array[j+1].mark.version);
            printf("Stuck5\n");
            continue;
          }
          if (pma->array[j+1].mark.operation == 0)
          {
            marker_t old = pma->array[j+1].mark;
            marker_t new = {.operation = 1, .key = 0, .val = old.key, .version = old.version + 1};
            if (CASM(&pma->array[j+1].mark, old, new))
            {
              uint64_t old_key = pma->array[j].key;
              uint64_t old_val = pma->array[j].val;
              if (DCAS(&pma->array[j].key, &pma->array[j].val, old_key, pma->array[j+1].key, old_val, pma->array[j+1].val))
              {
                pma->array[j].mark = pma->array[j+1].mark;
                pma->array[j].mark.operation = 0;
                pma->array[j].version = pma->array[j+1].mark.version;
                uint64_t old_key_from = pma->array[j+1].key;
                help = true;
                if (!CAS(&pma->array[j+1].key, old_key_from, 0ULL))
                {
                  // can fail here because can be helped
                }
                (pma->array[j+1].mark.operation = 0);
                (pma->array[j+1].version = pma->array[j+1].mark.version);
                break;
              }
              else
              {
                (pma->array[j+1].mark.operation = 0);
                (pma->array[j+1].version = pma->array[j+1].mark.version);
              }
            }
            else
            {
              printf("Stuck6\n");
            }
          }
          else
          {
            if (help)
            {
              pma->array[j+1].key = pma->array[j+1].mark.key;
              pma->array[j+1].val = pma->array[j+1].mark.val;
            }
          }
        }
        j++;
      }
      while (true)
      {
        bool help = false;
        if (pma->array[i - 1].version < pma->array[i - 1].mark.version)
        {
          // printf("Version: %d and Marker Version %d\n", pma->array[i-1].version ,pma->array[i-1].mark.version);
          printf("Stuck7\n");
          continue;
        }
        if (pma->array[i - 1].mark.operation == 0)
        {
          // perform CAS on the cell marker
          marker_t old = pma->array[i - 1].mark;
          marker_t new = {.operation = 1, .key = pma->array[i-1].key, .val = pma->array[i-1].val, .version = old.version + 1};
          if (CASM(&pma->array[i - 1].mark, old, new))
          {
            uint64_t old_key = pma->array[i].key;
            uint64_t old_val = pma->array[i].val;
            if (!DCAS(&pma->array[i].key, &pma->array[i].val, old_key, key, old_val, val))
            {
              pma->array[i - 1].mark.operation = 0;
              pma->array[i - 1].version = pma->array[i - 1].mark.version;
              continue;
            }
            pma->array[i - 1].mark.operation = 0;
            pma->array[i - 1].version = pma->array[i - 1].mark.version;
            break;
          }
          else
          {
            printf("Stuck8\n");
            continue;
          }
        }
      }
    }
  }
  pma->n++;
  rebalance(pma, i);
  // printf("New size of pma: %d\n", pma->m);
}

bool pma_delete(PMA pma, key_t key)
{
  int64_t i;
  if (pma_find(pma, key, &i))
  {
    pma_delete_at(pma, i);
    return (true);
  }
  return (false); /* key does not exist. */
}

void pma_delete_at(PMA pma, int64_t i)
{
  assert(i >= 0);
  assert(i < pma->m);
  keyval_clear(&(pma->array[i]));
  rebalance(pma, i);
}

void pma_get(PMA pma, int64_t i, keyval_t *keyval)
{
  assert(i >= 0);
  assert(i < pma->m);
  *keyval = pma->array[i];
}

void pma_print(PMA pma)
{
  for (int x = 0; x < pma->m; x++)
  {
    printf("Index: %d, key: %d, value: %d, version: %d\n", x, pma->array[x].key, pma->array[x].val, pma->array[x].version);
  }
}

/* Returns the size of the array. */
uint64_t pma_capacity(PMA p)
{
  return p->m;
}

/* Returns the number of elements in the array. */
uint64_t pma_count(PMA p)
{
  return p->n;
}

static void rebalance(PMA pma, int64_t i)
{
  // printf("index: %d\n",i);

  int64_t window_start, window_end;
  uint8_t height = 0;
  uint64_t occupancy = (keyval_empty(&(pma->array[i]))) ? 0 : 1;
  int64_t left_index = i - 1;
  int64_t right_index = i + 1;
  double density, t_height, p_height;
  do
  {
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
    density = (double)occupancy / (double)window_size;
    t_height = t_0 - (height * pma->delta_t);
    p_height = p_0 + (height * pma->delta_p);
    height++;
  } while ((density < p_height ||
            density >= t_height) &&
           height < pma->h);
  /* Found a window within threshold. */
  if (density >= p_height && density < t_height)
  {
    // printf("window_start: %d, window_end: %d\n", window_start, window_end);
    while (!(pack(pma, window_start, window_end, occupancy) &&
             spread(pma, window_start, window_end, occupancy)))
    {
      printf("FOUND WINDOW\n");
    };
  }
  else
  {
    while (!resize(pma))
      printf("THIS IS RESIZE\n");
  }
}

/* from is inclusive, to is exclusive. */
static bool pack(PMA pma, uint64_t from, uint64_t to, uint64_t n)
{
  // printf("FROM: %d & TO: %d , SUCCESS: %d\n", from, to, (from < to));
  assert(from < to);
  uint64_t read_index = from;
  uint64_t write_index = from;
  while (read_index < to)
  {
    if (!keyval_empty(&(pma->array[read_index])))
    {
      if (read_index > write_index)
      {
        bool help = false;
        // printf("version: %d and marker version: %d\n", pma->array[read_index].version, pma->array[read_index].mark.version);
        if (pma->array[read_index].version < pma->array[read_index].mark.version)
        {
          // printf("Version: %d and Marker Version %d\n", pma->array[read_index].version ,pma->array[read_index].mark.version);
          printf("Stuck10\n");
          return false;
        }
        if (pma->array[read_index].mark.operation == 0)
        {
          marker_t old = pma->array[read_index].mark;
          marker_t new = {.operation = 1, .key = 0, .val = old.key, .version = old.version + 1};
          if (CASM(&pma->array[read_index].mark, old, new))
          {
            uint64_t old_key = pma->array[write_index].key;
            uint64_t old_val = pma->array[write_index].val;
            if (DCAS(&pma->array[write_index].key, &pma->array[write_index].val, old_key, pma->array[read_index].key, old_val, pma->array[read_index].val))
            {
              pma->array[write_index].mark = pma->array[read_index].mark;
              pma->array[write_index].version = pma->array[read_index].mark.version;
              pma->array[write_index].mark.operation = 0;
              uint64_t old_key_from = pma->array[read_index].key;
              help = true;
              if (!CAS(&pma->array[read_index].key, old_key_from, 0ULL))
              {
                // can fail here because can be helped
              }
              (pma->array[read_index].mark.operation = 0);
              (pma->array[read_index].version = pma->array[read_index].mark.version);
            }
            else
            {
              (pma->array[read_index].mark.operation = 0);
              (pma->array[read_index].version = pma->array[read_index].mark.version);
              return false;
            }
          }
          else
          {
            printf("Stuck11\n");
            return false;
          }
        }
        else
        {
          if (help)
          {
            pma->array[read_index].key = pma->array[read_index].mark.key;
            pma->array[read_index].val = pma->array[read_index].mark.val;
          }
        }
      }
      write_index++;
    }
    read_index++;
  }
  return true;
  // assert (n == write_index - from);
}

/* from is inclusive, to is exclusive. */
static bool spread(PMA pma, uint64_t from, uint64_t to, uint64_t n)
{
  assert(from < to);
  uint64_t capacity = to - from;
  uint64_t frequency = (capacity << 8) / n; /* 8-bit fixed point arithmetic. */
  uint64_t read_index = from + n - 1;
  uint64_t write_index = (to << 8) - frequency;
  while ((write_index >> 8) > read_index)
  {
    bool help = false;
    if (pma->array[read_index].version < pma->array[read_index].mark.version)
    {
      // printf("Version: %d and Marker Version %d\n", pma->array[read_index].version ,pma->array[read_index].mark.version);
      printf("Stuck12\n");
      return false;
    }
    if (pma->array[read_index].mark.operation == 0)
    {
      marker_t old = pma->array[read_index].mark;
      marker_t new = {.operation = 1, .key = 0, .val = old.key, .version = old.version + 1};
      if (CASM(&pma->array[read_index].mark, old, new))
      {
        uint64_t old_key = pma->array[write_index >> 8].key;
        uint64_t old_val = pma->array[write_index >> 8].val;
        if (DCAS(&pma->array[write_index >> 8].key, &pma->array[write_index >> 8].val, old_key, pma->array[read_index].key, old_val, pma->array[read_index].val))
        {
          pma->array[write_index >> 8].mark = pma->array[read_index].mark;
          pma->array[write_index >> 8].mark.operation = 0;
          pma->array[write_index >> 8].version = pma->array[read_index].mark.version;
          uint64_t old_key_from = pma->array[read_index].key;
          help = true;
          if (!CAS(&pma->array[read_index].key, old_key_from, 0ULL))
          {
            // can fail here because can be helped
          }
          (pma->array[read_index].mark.operation = 0);
          (pma->array[read_index].version = pma->array[read_index].mark.version);
        }
        else
        {
          (pma->array[read_index].mark.operation = 0);
          (pma->array[read_index].version = pma->array[read_index].mark.version);
          return false;
        }
      }
      else
      {
        printf("Stuck13\n");
        return false;
      }
    }
    else
    {
      if (help)
      {
        pma->array[read_index].key = pma->array[read_index].mark.key;
      }
    }
    read_index--;
    write_index -= frequency;
  }
  return true;
}

static bool resize(PMA pma)
{
  if (!pack(pma, 0, pma->m, pma->n))
    return false;
  uint64_t old_m = pma->m;
  uint64_t new_m = compute_capacity(pma);
  pma->m = new_m;
  pma->h = floor_lg(pma->num_segments) + 1;
  pma->delta_t = (t_0 - t_h) / pma->h;
  pma->delta_p = (p_h - p_0) / pma->h;
  pma->array = (keyval_t *)realloc(pma->array, sizeof(keyval_t) * new_m);
  for (int x = old_m; x < new_m; x++)
  {
    pma->array[x].version = 0;
    marker_t mk1 = {.operation = 0, .key = 0, .val = 0, .version = 0};
    pma->array[x].mark = mk1;
  }
  for (uint64_t i = pma->n; i < pma->m; i++)
  {
    while (true)
    {
      bool help = false;
      if (pma->array[i-1].version < pma->array[i-1].mark.version)
      {
        // printf("Version: %d and Marker Version %d\n", pma->array[i-1].version ,pma->array[i-1].mark.version);
        printf("Stuck14\n");
        continue;
      }
      if (pma->array[i-1].mark.operation == 0)
      {
        // perform CAS on the cell marker
        marker_t old = pma->array[i-1].mark;
        marker_t new = {.operation = 1, .key = old.key, .val = old.val, .version = old.version + 1};
        if (CASM(&pma->array[i-1].mark, old, new))
        {
          uint64_t old_key = pma->array[i].key;
          uint64_t old_val = pma->array[i].val;
          if (!DCAS(&pma->array[i].key, &pma->array[i].val, old_key, 0ULL, old_val, 0ULL))
          {
            pma->array[i-1].mark.operation = 0;
            pma->array[i-1].version = pma->array[i-1].mark.version;
            continue;
          }
          pma->array[i-1].mark.operation = 0;
          pma->array[i-1].version = pma->array[i-1].mark.version;
          break;
        }
        else
        {
          printf("Stuck15\n");
          continue;
        }
      }
    }
  }

  if (!spread(pma, 0, pma->m, pma->n))
    return false;
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

void *threadFunc(void *args)
{
  int tid = *((int *)args);
  // printf("tid is %d \n", tid);
  for (uint64_t x = pow(10, tid); x <= pow(10, tid + 2); x++)
  {
    pma_insert(pma1, x, x);
  }
  pthread_exit(NULL);
}
int main(int argc, char *argv[])
{
  pma1 = pma_create();
  pthread_t tid[7];
  int msec = 0, trigger = 10; /* 10ms */
  clock_t before = clock();
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
  // printf("after pthread_join \n");
  clock_t difference = clock() - before;
  msec = difference * 1000 / CLOCKS_PER_SEC;
  printf("Time taken: %d seconds %d milliseconds\n",
         msec / 1000, msec % 1000);
  // pma_print(pma1);
  printf("Elements: %d\n", pma1->n);
  printf("Capacity: %d", pma1->m);
  pthread_exit(NULL);
}
