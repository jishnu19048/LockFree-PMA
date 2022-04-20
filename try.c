#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <pthread.h>

typedef struct
{
  int operation;
  uint64_t version;
  int key;
  int val;
} marker_t;

int main() {
    marker_t x = {.key=4, .operation=5, .val=6, .version=7};
    marker_t y = x;
    y.key = 7;
    printf("%d %d %d %d\n", y.key, y.val, y.operation, y.version);
    printf("%d %d %d %d\n", x.key, x.val, x.operation, x.version);
    return 0;
}