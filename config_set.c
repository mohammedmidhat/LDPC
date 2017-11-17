#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>


#define Q		5


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
	}/*
	if(rest){
		for(i = 0; i <= Q-rest; i++){
			conv_sequence[seq_sz_total - seq_sz] = i%Q;
			find_all_config(seq_sz_total, seq_sz-1, rest+i%Q, not_row);
		}
	} else{
		for(i = 0; i < Q; i++){
			conv_sequence[seq_sz_total - seq_sz] = i%Q;
			find_all_config(seq_sz_total, seq_sz-1, rest+i%Q, not_row);
		}
	}*/
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


void main(void){
	int i, j;
	int rmax = 3;
	conv_sequence = malloc((rmax-1)*sizeof(int));

	fQa[0][0] = 0.45;
	fQa[0][1] = 0.25;
	fQa[0][2] = 0.15;
	fQa[0][3] = 0.05;
	fQa[0][4] = 0.1;
	fQa[1][0] = 0.23;
	fQa[1][1] = 0.64;
	fQa[1][2] = 0.1;
	fQa[1][3] = 0.02;
	fQa[1][4] = 0.01;
	fQa[2][0] = 0.22;
	fQa[2][1] = 0.29;
	fQa[2][2] = 0.46;
	fQa[2][3] = 0.03;
	fQa[2][4] = 0;
	fQa[3][0] = 0.22;
	fQa[3][1] = 0.29;
	fQa[3][2] = 0.46;
	fQa[3][3] = 0.03;
	fQa[3][4] = 0;

	CNP(4,4);

	for(i = 0; i < rmax; i++){
		printf("fRa[%d] = ", i);
		for(j = 0; j < Q; j++){
			printf("%.2f ", fRa[i][j]);
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
