#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <errno.h>
#include <time.h>
 

#define Q	13
#define Log2Q	4
#define half_Q_binary_ext  1<<(Log2Q-1)

float ntt_temp[half_Q_binary_ext<<1];

void ntt(float p[Q])
{
  memcpy(ntt_temp, p, Q*sizeof(float));

  int b, factor = 1, rest;
  for (b = 0; b < Log2Q; b++) {
    for (rest = 0; rest < 8; rest++) {
      int restH = rest >> b;
      int restL = rest & (factor-1);
      int rest0 = (restH << (b+1)) + restL;
      int rest1 = rest0 + factor;
      float prest0 = ntt_temp[rest0];

      ntt_temp[rest0] += ntt_temp[rest1];
      ntt_temp[rest1] = prest0 - ntt_temp[rest1];
    }
    factor += factor;
  }

  memcpy(p, ntt_temp, Q*sizeof(float));
}

int main(void){
	int i;
	float p[Q] = {0.4,3.5,6.7,0.8,9.8,7.6,5.4,3.35,6.7,8.9,8.7,0.6,0.5};
	memset(ntt_temp+Q, 0, (16-Q)*sizeof(float));

	for(i = 0; i < Q; i++){
		printf("%.2f ", p[i]);
	}
	printf("\n");

	ntt(p);

	for(i = 0; i < Q; i++){
		printf("%.2f ", p[i]);
	}
	printf("\n");

	ntt(p);

	for(i = 0; i < Q; i++){
		printf("%.2f ", p[i]/16.0);
	}
	printf("\n");

	ntt(p);

	for(i = 0; i < Q; i++){
		printf("%.2f ", p[i]/16.0);
	}
	printf("\n");

	return 0;
}
