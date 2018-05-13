#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <errno.h>
#include <time.h>
#include "tajj.h"

#define INT   8
#define FRAC  8
#define FRAC_LEVELS 256

void display_two_comp(int val){
  int i;
  for(i = INT+FRAC-1; i >= 0; i--){
    printf("%d", (val*FRAC_LEVELS >> i)&1);
  }
  printf("\n");
}

void main(void){
  double y = 42.793;
  int i;
  for(i = INT+FRAC-1; i >= 0; i--){
    printf("%d", ((int)(y*FRAC_LEVELS) >> i)&1);
  }
}
