#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>


#define Q		4
#define Log2Q	2


int *conv_sequence;
float p_symbol_alloc;
float *p_symbol = &p_symbol_alloc;
float fQa[4][Q];
float fRa[4][Q];


void find_all_config(int seq_sz_total, int seq_sz, int rest, int not_row){
  //printf("%d, %d, %d, %d\n", seq_sz_total, seq_sz, rest, not_row);
  int i, j;
  int acc;
  
  if(seq_sz == 1){
    int seq_iter = 0;
    float p_seq = 1.0;
    
    conv_sequence[seq_sz_total - seq_sz] = ((Q - rest)%Q+Q)%Q;
    printf("seq = ");
    for(i = 0; i < seq_sz_total+1; i++){
      if(i != not_row){
        printf("%d ", conv_sequence[seq_iter]);
        p_seq *= fQa[i][conv_sequence[seq_iter]];
        seq_iter++;
      }
    }
    *p_symbol += p_seq;
    printf("\n");
    return;
  }
  for(i = 0; i < Q; i++){
    conv_sequence[seq_sz_total - seq_sz] = i;
    find_all_config(seq_sz_total, seq_sz-1, rest+i, not_row);
  }
}

void CNP(int row_weight, int syndrome){
	int i, a;
	
	for(i = 0; i < row_weight; i++){
		printf("new config\n");
		for (a = 0; a < Q; a++){
			printf("config[%d]\n", a);
			*p_symbol = 0;
			find_all_config(row_weight-1, row_weight-1, (((a-syndrome)%Q)+Q)%Q, i);
			fRa[i][a] = *p_symbol;
		}
	}
}

void ntt(float p[Q])
{
  int b, factor = 1, rest;
  for (b = 0; b < Log2Q; b++) {
    for (rest = 0; rest < (Log2Q<<1)/2; rest++) {
      int restH = rest >> b;
      int restL = rest & (factor-1);
      int rest0 = (restH << (b+1)) + restL;
      int rest1 = rest0 + factor;
      float prest0 = p[rest0];

      p[rest0] += p[rest1];
      p[rest1] = prest0 - p[rest1];
    }
    factor += factor;
  }
}


void main(void){
	int i, j, k;
	int rmax = 3;
	conv_sequence = malloc((rmax-1)*sizeof(int));

	fQa[0][0] = 0.8;
	fQa[0][1] = 0.05;
	fQa[0][2] = 0.1;
	fQa[0][3] = 0.05;
	//fQa[0][3] = 0.0;

	fQa[1][0] = 0.03;
	fQa[1][1] = 0.05;
	fQa[1][2] = 0.9;
	fQa[1][3] = 0.02;
	//fQa[1][3] = 0.0;

	fQa[2][0] = 0.05;
	fQa[2][1] = 0.89;
	fQa[2][2] = 0.03;
	fQa[2][3] = 0.03;
	//fQa[2][3] = 0.0;

	/*fQa[3][0] = 0.02;
	fQa[3][1] = 0.02;
	fQa[3][2] = 0.11;
	fQa[3][3] = 0.85;*/

	for(i = 0; i < rmax; i++){
		printf("fqa[%d] = ", i);
		for(j = 0; j < (Log2Q<<1); j++){
			printf("%.2f ", fQa[i][j]);
		}
		printf("\n");
	}

	CNP(rmax,0);

	printf("convolution:\n");
	for(i = 0; i < rmax; i++){
		printf("fRa[%d] = ", i);
		for(j = 0; j < (Log2Q<<1); j++){
			printf("%.2f ", fRa[i][j]);
		}
		printf("\n");
	}

	printf("fft:\n");
	for(i = 0; i < rmax; i++){
		ntt(fQa[i]);
	}

	for(i = 0; i < rmax; i++){
		printf("fQa[%d] = ", i);
		for(j = 0; j < (Log2Q<<1); j++){
			printf("%.2f ", fQa[i][j]);
		}
		printf("\n");
	}
	printf("\n");

	for(i = 0; i < rmax; i++){
		for(j = 0; j < (Log2Q<<1); j++){
			fRa[i][j] = 1.0;
		}
		for(j = 0; j < rmax; j++){
			if(i != j){
				for(k = 0; k < (Log2Q<<1); k++){
					fRa[i][k] *= fQa[j][k];
				}
			}
		}
	}

	for(i = 0; i < rmax; i++){
		printf("fRa[%d] = ", i);
		for(j = 0; j < (Log2Q<<1); j++){
			printf("%.2f ", fRa[i][j]);
		}
		printf("\n");
	}
	printf("\n");

	for(i = 0; i < rmax; i++){
		ntt(fRa[i]);
	}

	for(i = 0; i < rmax; i++){
		printf("fRa[%d] = ", i);
		for(j = 0; j < (Log2Q<<1); j++){
			printf("%.2f ", fRa[i][j]/4.0);
		}
		printf("\n");
	}
/*
	FILE *fp = fopen("sn.bin", "rb");
	char buf[3];
	fgets(buf, sizeof(buf), fp);
	printf("hi\n");
	for(int i = 0; i < 3; i++){
		printf("%u ", buf[i]);
		printf("hi\n");
	}*/
}
