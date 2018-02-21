#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <errno.h>
#include <time.h>
#include "tajj.h"

#define INT   4
#define FRAC    4
#define INT_LEVELS  pow(2,INT-1)
#define FRAC_LEVELS   pow(2,FRAC)

double float_to_fix(double val){
  double intermed = round(val*FRAC_LEVELS);
  double result = intermed/FRAC_LEVELS;
  if(result > INT_LEVELS-1){
    return INT_LEVELS-1;
  } else if(result < -INT_LEVELS){
    return -INT_LEVELS;
  } else{
    return result;
  }
}

void main(void){
	int i;
	double x[5] = {8, -8, -8.125, 2.35, -4.31};

	for(i = 0; i < 5; i++){
		printf("%.3f ", float_to_fix(x[i]));
	}
}
