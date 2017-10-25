#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>

#define alphabet_sz		7
#define seq_total		4

int y[seq_total];

void find_all_config(int seq_sz, int rest){
	int i, j;
	int acc;
	
	if(seq_sz == 1){
		y[seq_total - seq_sz] = (alphabet_sz - rest)%alphabet_sz;
		printf("seq = ");
		for(i = 0; i < seq_total; i++){
			printf("%d ", y[i]);
		}
		printf("\n");
		return;
	}
	if(rest){
		for(i = 0; i <= alphabet_sz-rest; i++){
			acc = 0;

			y[seq_total - seq_sz] = i%alphabet_sz;
			for(j = 0; j < seq_total - seq_sz + 1; j++){
				acc+= y[j];
			}
		
			find_all_config(seq_sz-1, acc);
		}
	} else{
		for(i = 0; i < alphabet_sz; i++){
			acc = 0;

			y[seq_total - seq_sz] = i%;
			for(j = 0; j < seq_total - seq_sz + 1; j++){
				acc+= y[j];
			}
		
			find_all_config(seq_sz-1, acc);
		}
	}
}

void main(void){
	find_all_config(seq_total,0);
}