#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <errno.h>
#include <time.h>


#define Q 13
#define page_size 18336
#define QLC   4
#define CW_per_page   1


int grey_code_inv[16] = {6,5,7,14,9,12,8,13,3,4,2,15,10,11,1,0};
int n;
// Lower page data stored first
char data_read[QLC*page_size];
char data_written[QLC*page_size];


void get_symbols_in_byte(int x[], int byte_ind, char lp_byte, char mp_byte, char up_byte, char tp_byte){
  int i, symbol, lp_bit, mp_bit, up_bit, tp_bit;
  for(i = 0; i < 8; i++){
    lp_bit = (lp_byte >> (8-1-i)) & 1;
    mp_bit = (mp_byte >> (8-1-i)) & 1;
    up_bit = (up_byte >> (8-1-i)) & 1;
    tp_bit = (tp_byte >> (8-1-i)) & 1;

    symbol = grey_code_inv[(tp_bit << 3) + (up_bit << 2) + (mp_bit << 1) + (lp_bit)];
    x[8*byte_ind+i] = symbol;
  }
}


void main(void){
    FILE *written_data = fopen("written_data.bin", "rb");

    FILE *sym = fopen("sym.txt", "w");
    
    for(i = 0; i < num_trials/CW_per_page; i++){
      fread(data_written, 1, QLC*page_size, written_data);

      for(j = 0; j < CW_per_page; j++){
        for(k = 0; k < block_sz_bytes; k++){
          lp_byte = data_written[j*block_sz_bytes+k];
          mp_byte = data_written[page_size+j*block_sz_bytes+k];
          up_byte = data_written[2*page_size+j*block_sz_bytes+k];
          tp_byte = data_written[3*page_size+j*block_sz_bytes+k];
          get_symbols_in_byte(x, k, lp_byte, mp_byte, up_byte, tp_byte);
        }
        channel(x, y, p_rec_given_sent, logfna, num_reads);
        for(k = 0; k < n; k++){
          fprintf(sym, "%d ", y[k]);
        }


        /*enc(x, s);
        assign_llr(y, p_rec_given_sent, logfna, num_reads);
        dec_result = dec(logfna, s, iteration);
    
        if(dec_result){
          errors[0]++;
        } else {
          if(HamDist(tmp_x, x, n) != 0) errors[1]++;
        }*/
      }
    }

    //fclose(read_data);
    fclose(written_data);

    fclose(sym);
}
