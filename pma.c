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

/*
 * Returns the 1-based index of the last (i.e., most significant) bit set in x.
 */
uint64_t last_bit_set (uint64_t x) {
  // assert (x > 0);
  return (sizeof (uint64_t) * 8 - __builtin_clzll (x));
}

uint64_t floor_lg (uint64_t x) {
  return (last_bit_set (x) - 1);
}

uint64_t ceil_lg (uint64_t x) {
  return (last_bit_set (x - 1));
}

/* Returns the largest power of 2 not greater than x
 * (i.e., $2^{\lfloor \lg x \rfloor}$).
 */
uint64_t hyperfloor (uint64_t x) {
  return (1 << floor_lg (x));
}

/* Returns the smallest power of 2 not less than x
 * (i.e., $2^{\lceil \lg x \rceil}$).
 */
uint64_t hyperceil (uint64_t x) {
  return (1 << ceil_lg (x));
}

uint64_t ceil_div (uint64_t x, uint64_t y) {
  assert (x > 0);
  return (1 + ((x - 1) / y));
}

typedef uint64_t key_t;
typedef uint64_t val_t;

/* compare and swap implementation */
bool CAS( int *p, int old, int new ) {
  if( *p == old ) {
    *p = new;
    return true;
  } else {
    return false;
  }
}
typedef struct{
  int operation;
  uint64_t version;
  key_t key;
  val_t val;
} marker_t;
typedef struct {
  key_t key;
  val_t val;
  marker_t mark;
  uint64_t version;
} keyval_t;

/* Returns true if keyval is empty and false otherwise. */
bool keyval_empty (const keyval_t *keyval) {
  return (keyval->key == 0ULL);
}

/* Sets keyval to be empty. */
void keyval_clear (keyval_t *keyval) {
  keyval->key = 0ULL;
  keyval->val = 0ULL;
}

/* Reserve 8 bits to allow for fixed point arithmetic. */
#define MAX_SIZE ((1ULL << 56) - 1ULL)

/* Height-based (as opposed to depth-based) thresholds. */
/* Upper density thresholds. */
static const double t_h = 0.75;  /* root. */
static const double t_0 = 1.00;  /* leaves. */
/* Lower density thresholds. */
static const double p_h = 0.50;  /* root. */
static const double p_0 = 0.25;  /* leaves. */

static const uint8_t max_sparseness = 1 / p_0;
static const uint8_t largest_empty_segment = 1 * max_sparseness;

struct _pma;
typedef struct _pma pma_t, *PMA;

PMA pma_create (void);
PMA pma_from_array (keyval_t *array, uint64_t n);
void pma_destroy (PMA *pma);
bool pma_find (PMA pma, key_t key, int64_t *index);
bool pma_insert (PMA pma, key_t key, val_t val);
void pma_insert_after (PMA pma, int64_t i, key_t key, val_t val);
bool pma_delete (PMA pma, key_t key);
void pma_delete_at (PMA pma, int64_t i);
void pma_get (PMA pma, int64_t i, keyval_t *keyval);
uint64_t pma_capacity (PMA p);
uint64_t pma_count (PMA p);

/* TODO: For testing purposes only. */
uint8_t pma_segment_size (PMA pma);


struct _pma {
  uint64_t n;  /* Number of elements. */
  uint64_t m;  /* Size of the array (capacity). */
  uint8_t s;  /* Size of the segments. */
  uint64_t num_segments;
  uint8_t h;  /* Height of the tree. */
  double delta_t;  /* Delta for upper density threshold. */
  double delta_p;  /* Delta for lower density threshold. */
  keyval_t *array;
};

/* TODO: For testing purposes only. */
uint8_t pma_segment_size (PMA pma) {
  return (pma->s);
}

static void rebalance (PMA pma, int64_t i);
static void pack (PMA pma, uint64_t from, uint64_t to, uint64_t n);
static void spread (PMA pma, uint64_t from, uint64_t to, uint64_t n);
static void resize (PMA pma);
static void compute_capacity (PMA pma);

PMA pma_create () {
  PMA pma = (PMA)malloc (sizeof (pma_t));
  pma->n = 0;
  pma->s = largest_empty_segment;
  pma->m = (1ULL << largest_empty_segment);
  pma->num_segments = pma->m / pma->s;
  pma->h = floor_lg (pma->num_segments) + 1;
  pma->delta_t = (t_0 - t_h) / pma->h;
  pma->delta_p = (p_h - p_0) / pma->h;
  pma->array = (keyval_t *)malloc (sizeof (keyval_t) * pma->m);
  return (pma);
}

PMA pma_from_array (keyval_t *array, uint64_t n) {
  assert (n > 0);
  PMA pma = (PMA)malloc (sizeof (pma_t));
  pma->n = n;
  compute_capacity (pma);
  pma->h = floor_lg (pma->num_segments) + 1;
  pma->delta_t = (t_0 - t_h) / pma->h;
  pma->delta_p = (p_h - p_0) / pma->h;
  pma->array = (keyval_t *)malloc (sizeof (keyval_t) * pma->m);
  for (uint64_t i = 0; i < pma->m; i++) {
    if (i < pma->n) {
      pma->array [i].key = array [i].key;
      pma->array [i].val = array [i].val;
    } else {
      keyval_clear (&(pma->array [i]));
    }
  }
  spread (pma, 0, pma->m, pma->n);
  return (pma);
}

void pma_destroy (PMA *pma) {
  free ((*pma)->array);
  (*pma)->array = NULL;
  free (*pma);
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
bool pma_find (PMA pma, key_t key, int64_t *index) {
  int64_t from = 0;
  int64_t to = pma->m - 1;
  while (from <= to) {
    int64_t mid = from + (to - from) / 2;
    int64_t i = mid;
    /* Start scanning left until we find a non-empty slot or we reach past the
     * beginning of the subarray. */
    while (i >= from &&
           keyval_empty (&(pma->array [i])))
      i--;
    if (i < from) {  /* Everything between from and mid (inclusive) is empty. */
      from = mid + 1;
    } else {
      if (pma->array [i].key == key) {
        *index = i;
        return (true);
      }
      else if (pma->array [i].key < key) {
        from = mid + 1;
      } else {  /* pma->array [i].key > key */
        to = i - 1;
      }
    }
  }
  /* Didn't find `key'. `to' should hold its predecessor (unless it's empty). */
  *index = to;
  while (*index >= 0 && keyval_empty (&(pma->array [*index])))
    (*index)--;
  return (false);
}

bool pma_insert (PMA pma, key_t key, val_t val) {
  int64_t i;
  if (!pma_find (pma, key, &i)) {  /* We do not allow duplicates.*/
    pma_insert_after (pma, i, key, val);
    return (true);
  }
  return (false);
}

void pma_insert_after (PMA pma, int64_t i, key_t key, val_t val) {
  assert (i >= -1);
  assert (i < (int64_t)pma->m);
  assert (i >= 0 && !keyval_empty (&(pma->array [i])) || i >= -1);
  /* Find an empty space to the right of i. There should be one close by. */
  int64_t j = i + 1;
  while (j < pma->m && !keyval_empty (&(pma->array [j])))
    j++;
  if (j < pma->m) {  /* Found one. */
    while (j > i + 1) {  /* Push elements to make space for the new element. */
      pma->array [j].key = pma->array [j - 1].key;
      pma->array [j].val = pma->array [j - 1].val;
      j--;
    }
    pma->array [i + 1].key = key;
    pma->array [i + 1].val = val;
    i = i + 1;  /* Update i to point to where the new element is. */
  } else {  /* No empty space to the right of i. Try left. */
    j = i - 1;
    while (j >= 0 && !keyval_empty (&(pma->array [j])))
      j--;
    if (j >= 0) {  /* Found one. */
      while (j < i) {  /* Push elements to make space for the new element. */
        pma->array [j].key = pma->array [j + 1].key;
        pma->array [j].val = pma->array [j + 1].val;
        j++;
      }
      pma->array [i].key = key;
      pma->array [i].val = val;
    }
  }
  pma->n++;
  rebalance (pma, i);
}

bool pma_delete (PMA pma, key_t key) {
  int64_t i;
  if (pma_find (pma, key, &i)) {
    pma_delete_at (pma, i);
    return (true);
  }
  return (false);  /* key does not exist. */
}

void pma_delete_at (PMA pma, int64_t i) {
  assert (i >= 0);
  assert (i < pma->m);
  keyval_clear (&(pma->array [i]));
  rebalance (pma, i);
}

void pma_get (PMA pma, int64_t i, keyval_t *keyval) {
  assert (i >= 0);
  assert (i < pma->m);
  *keyval = pma->array [i];
}

/* Returns the size of the array. */
uint64_t pma_capacity (PMA p) {
  return p->m;
}

/* Returns the number of elements in the array. */
uint64_t pma_count (PMA p) {
  return p->n;
}

static void rebalance (PMA pma, int64_t i) {
  int64_t window_start, window_end;
  uint8_t height = 0;
  uint64_t occupancy = (keyval_empty (&(pma->array [i]))) ? 0 : 1;
  int64_t left_index = i - 1;
  int64_t right_index = i + 1;
  double density, t_height, p_height;
  do {
    uint64_t window_size = pma->s * (1 << height);
    uint64_t window = i / window_size;
    window_start = window * window_size;
    window_end = window_start + window_size;
    while (left_index >= window_start) {
      if (!keyval_empty (&(pma->array [left_index])))
        occupancy++;
      left_index--;
    }
    while (right_index < window_end) {
      if (!keyval_empty (&(pma->array [right_index])))
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
  if (density >= p_height && density < t_height) {
    pack (pma, window_start, window_end, occupancy);
    spread (pma, window_start, window_end, occupancy);
  } else {
    resize (pma);
  }
}

/* from is inclusive, to is exclusive. */
static void pack (PMA pma, uint64_t from, uint64_t to, uint64_t n) {
  assert (from < to);
  uint64_t read_index = from;
  uint64_t write_index = from;
  while (read_index < to) {
    if (!keyval_empty (&(pma->array [read_index]))) {
      if (read_index > write_index) {
        pma->array [write_index].key = pma->array [read_index].key;
        pma->array [write_index].val = pma->array [read_index].val;
        keyval_clear (&(pma->array [read_index]));
      }
      write_index++;
    }
    read_index++;
  }
  // assert (n == write_index - from);
}

/* from is inclusive, to is exclusive. */
static void spread (PMA pma, uint64_t from, uint64_t to, uint64_t n) {
  assert (from < to);
  uint64_t capacity = to - from;
  uint64_t frequency = (capacity << 8) / n;  /* 8-bit fixed point arithmetic. */
  uint64_t read_index = from + n - 1;
  uint64_t write_index = (to << 8) - frequency;
  while ((write_index >> 8) > read_index) {
    pma->array [write_index >> 8].key = pma->array [read_index].key;
    pma->array [write_index >> 8].val = pma->array [read_index].val;
    keyval_clear (&(pma->array [read_index]));
    read_index--;
    write_index -= frequency;
  }
}

static void resize (PMA pma) {
  pack (pma, 0, pma->m, pma->n);
  compute_capacity (pma);
  pma->h = floor_lg (pma->num_segments) + 1;
  pma->delta_t = (t_0 - t_h) / pma->h;
  pma->delta_p = (p_h - p_0) / pma->h;
  pma->array = (keyval_t *)realloc (pma->array, sizeof (keyval_t) * pma->m);
  for (uint64_t i = pma->n; i < pma->m; i++)
    keyval_clear (&(pma->array [i]));
  spread (pma, 0, pma->m, pma->n);
}

static void compute_capacity (PMA pma) {
  pma->s = ceil_lg (pma->n);  /* Ideal segment size. */
  pma->num_segments = ceil_div (pma->n, pma->s);  /* Ideal number of segments. */
  /* The number of segments has to be a power of 2. */
  pma->num_segments = hyperceil (pma->num_segments);
  /* Update the segment size accordingly. */
  pma->s = ceil_div (pma->n, pma->num_segments);
  pma->m = pma->s * pma->num_segments;
  /* Scale up as much as possible. */
  pma->m = max_sparseness * pma->m;
  pma->s = max_sparseness * pma->s;
  assert (pma->m <= MAX_SIZE);
  assert (pma->m > pma->n);
}

int main(int argc, char *argv[]) {
  PMA pma1 = pma_create();
  int msec = 0, trigger = 10; /* 10ms */
  clock_t before = clock();
  for(int i=0;i<30;i++){
    pma_insert(pma1, i+1, i+1);
    clock_t difference = clock() - before;
    msec = difference * 1000 / CLOCKS_PER_SEC;
  }
  printf("Time taken: %d seconds %d milliseconds\n",
  msec/1000, msec%1000);
  printf("Elements: %d\n", pma1->n);
  printf("Capacity: %d", pma1->m);
}
