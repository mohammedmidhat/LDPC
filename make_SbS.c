#include <math.h>
#include <stdlib.h>
#include <stdio.h>


#define page_size	18336
#define QLC		4
#define Q	13
#define NUM_STR		4*64


int grey_code_inv[16] = {6,5,7,14,9,12,8,13,3,4,2,15,10,11,1,0};

double conf_mat[Q][Q];
// Lower page data stored first
char data_read[QLC*page_size];
char data_written[QLC*page_size];


void main(void){
	int i, j, k, symbol_written;
	int lp_bit_written, mp_bit_written, up_bit_written, tp_bit_written;
	int lp_bit_read, mp_bit_read, up_bit_read, tp_bit_read;
	char lp_byte_written, mp_byte_written, up_byte_written, tp_byte_written;
	char lp_byte_read, mp_byte_read, up_byte_read, tp_byte_read;
	double normal_factor = NUM_STR*page_size/13.0;

	// Initialize the zero values in the conf_mat
	for(i = 0; i < Q; i++){
		for(j = 0; j < Q; j++){
			conf_mat[i][j] = 0;
		}
	}

	FILE *read_data = fopen("read_data.bin", "rb");
    FILE *written_data = fopen("written_data.bin", "rb");

    for(i = 0; i < NUM_STR; i++){
    	fread(data_read, 1, QLC*page_size, read_data);
      	fread(data_written, 1, QLC*page_size, written_data);

      	for(j = 0; j < page_size; j++){
      		lp_byte_written = data_written[j];
          	mp_byte_written = data_written[page_size+j];
          	up_byte_written = data_written[2*page_size+j];
          	tp_byte_written = data_written[3*page_size+j];

          	lp_byte_read = data_read[j];
          	mp_byte_read = data_read[page_size+j];
          	up_byte_read = data_read[2*page_size+j];
          	tp_byte_read = data_read[3*page_size+j];

          	// Fetching the symbols from 4 bytes, 1 per QLC page
          	for(k = 0; k < 8; k++){
          		lp_bit_written = (lp_byte_written >> (8-1-k)) & 1;
    			mp_bit_written = (mp_byte_written >> (8-1-k)) & 1;
    			up_bit_written = (up_byte_written >> (8-1-k)) & 1;
    			tp_bit_written = (tp_byte_written >> (8-1-k)) & 1;

    			symbol_written = grey_code_inv[(tp_bit_written << 3) + (up_bit_written << 2) + (mp_bit_written << 1) + (lp_bit_written)];

    			lp_bit_read = (lp_byte_read >> (8-1-k)) & 1;
    			mp_bit_read = (mp_byte_read >> (8-1-k)) & 1;
    			up_bit_read = (up_byte_read >> (8-1-k)) & 1;
    			tp_bit_read = (tp_byte_read >> (8-1-k)) & 1;

    			symbol_read = grey_code_inv[(tp_bit_read << 3) + (up_bit_read << 2) + (mp_bit_read << 1) + (lp_bit_read)];

    			conf_mat[symbol_read][symbol_written]++;
          	}
      	}
    }

    // Normalizing conf_mat
	for(i = 0; i < Q; i++){
		for(j = 0; j < Q; j++){
			conf_mat[i][j] /= normal_factor;
			printf("%.3f ", conf_mat[i][j]);
		}
		printf("\n");
	}
}
