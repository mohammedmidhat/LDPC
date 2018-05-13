// Binary LDPC decoder, called in from a MATLAB script
// Uses Min-Sum algorithm
// With fixed-point math

// 3/30/2018
// Modified by Mohammed Al Ai Baky from Seishi Takamura's LDPC decoder implemention


#include "mex.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <limits.h>


// Each message is represented as 2's complement 
// with (INT) bits for the integer part (including the sign)
// and (FRAC) bits for the fractional part
#define INT   8
#define FRAC    8
int INT_LEVELS;
int FRAC_LEVELS;


int n, m;
int rmax, cmax;
int *row_weight, *col_weight;
int **row_col;

int **col_row;   // used in initdec() and for the Verilog decoder debugging in dec()
int *var_neighbors_displayed;   // allocated in initdec()
                               // the number messages to var neighbors already displayed
                               // used for the Verilog decoder debugging in dec()

// For Verilog decoder debugging
FILE *Verilog_sim_f;
int circ_size = 128;   // change for different codes

// Emulates fixed-point values
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

// check if there is at least one element in arr belongs to [low, high]
// low and high included
// returns 1 if this condition is satisfied
// return 0 otherwise
// len: is the array length
int array_element_in_range(int *arr, int len, int low, int high){
  int i;

  for(i = 0; i < len; i++){
    if((arr[i] >= low) && (arr[i] <= high)){
      return 1;
    }
  }

  return 0;
}

// For Verilog decoder debugging
void display_two_comp(double val){
  printf("%04x", (int)(val*FRAC_LEVELS));
  printf("\n");
}

double atanh2(double x)
{
  return log((1.0 + x) / (1.0 - x));  // returns 2*atanh(x)
}
double logtanh2(double x)
{
  return log(tanh(fabs(x*0.5)));  // returns log tanh |x|
}

int HamDist(int *x, int *y, int len)
{
  int i, sum = 0;
  for (i = 0; i < len; i++) {
    if (*x++ != *y++) sum++;
  }
  return sum;
}

int bsc(int x[], int y[], double p, double q0[])
{
  int i, num = 0, modified = 0;
  int *err = malloc(sizeof(int) * n);
  memset(err, 0, sizeof(int) * n);
  modified = n * p + 0.5;
  p = modified / (double)n; // correct error probability
  //printf("m/n=%g, ", (double)m/n);
  //printf("BSC channel entropy(rate) = %g (bits)\n",
  //       (-p*log(p)-(1-p)*log(1-p)) / log(2.0));
  while (num < modified) {
    i = rand() % n;
    if (err[i] == 1) continue;
    err[i] = 1;
    num++;
  }
  for (i = 0; i < n; i++) {
    y[i] = x[i] ^ err[i];
  }
  free(err);

  for (i = 0; i < n; i++) {
    double d = (1 - 2 * y[i]) * log((1.0 - p) / p);
    q0[i] = float_to_fix(d);
  }
  return modified;
}

void enc(int y[], int s[])
{
  int i, j;
  for (j = 0; j < m; j++) {
    register int k = 0;
    for (i = 0; i < row_weight[j]; i++)
      k ^= y[row_col[j][i]];

    s[j] = k;
  }
}

int **malloc2Dint(int a, int b) // allocates array[a][b]
{
  int i;
  int **pp = malloc(sizeof(int *) * a);
  int *p = malloc(sizeof(int) * a * b);
  if (pp == NULL || p == NULL) exit(-1);
  for (i = 0; i < a; i++) {
    pp[i] = p + b*i;
  }
  return pp;
}

int ***malloc2Dintp(int a, int b) // allocates array[a][b]
{
  int i;
  int ***pp = malloc(sizeof(int **) * a);
  int **p = malloc(sizeof(int*) * a * b);
  if (pp == NULL || p == NULL) exit(-1);
  for (i = 0; i < a; i++) {
    pp[i] = p + b*i;
  }
  return pp;
}

double **malloc2Ddouble(int a, int b) // allocates array[a][b]
{
  int i;
  double **pp = malloc(sizeof(double *) * a);
  double *p = malloc(sizeof(double) * a * b);
  if (pp == NULL || p == NULL) exit(-1);
  for (i = 0; i < a; i++) {
    pp[i] = p + b*i;
  }
  return pp;
}

double ***malloc2Ddoublep(int a, int b) // allocates array[a][b]
{
  int i;
  double ***pp = malloc(sizeof(double **) * a);
  double **p = malloc(sizeof(double*) * a * b);
  if (pp == NULL || p == NULL) exit(-1);
  for (i = 0; i < a; i++) {
    pp[i] = p + b*i;
  }
  return pp;
}


double **qin, ***qin_row;
double **LogTanhtin, ***LogTanhtin_row;
int **Sgntin, ***Sgntin_row;
int *tmp_bit;
int *tmp_s;


int dec(double q0[], int s[], int loop_max, int x[])
{
  int i, j, k, l, loop;
  
  memset(*qin, 0, n * cmax * sizeof(double));

  // Check if received data satisfied the code
  for (i = 0; i < n; i++) {
    tmp_bit[i] = (q0[i] < 0) ? 1 : 0;
  }
  enc(tmp_bit, tmp_s);
  i = HamDist(s, tmp_s, m);
  if (i == 0){
    fprintf(Verilog_sim_f, "already cw\n");     // For Verilog decoder debugging
    return 0;
  }


  for (loop = 0; loop < loop_max; loop++) {
    memset(var_neighbors_displayed, 0, n * sizeof(int));

    for (i = 0; i < n; i++) {
      double sum = q0[i];
      for (j = 0; j < col_weight[i]; j++){
        sum += qin[i][j];
      }

      for (j = 0; j < col_weight[i]; j++) {
        double qout = sum - qin[i][j];
        
        if (qout < 0) {
          *LogTanhtin_row[i][j] = float_to_fix(-qout);
          *Sgntin_row[i][j] = 1;
        } else {
          *LogTanhtin_row[i][j] = float_to_fix(qout);
          *Sgntin_row[i][j] = 0;
        }
      }
      // if(i == 0){
      //   printf("q0 = %08x\n", (int)(q0[i]*FRAC_LEVELS));
      //   for (j = 0; j < col_weight[0]; j++){
      //     printf("%08x ", (int)(qin[0][j]*FRAC_LEVELS));
      //   }
      //   printf("\n");
      //   printf("sum = %08x\n", (int)(sum*FRAC_LEVELS));
      //   for (j = 0; j < col_weight[0]; j++){
      //     printf("%08x ", (int)((sum - qin[0][j])*FRAC_LEVELS));
      //   }
      //   printf("\n\n");
      // }
    }

    // For Verilog decoder debugging (for QC-code only)
    for(j = 0; j < circ_size; j++){
      for(i = 0; i < m/circ_size; i++){
      //for(i = 0; i < cmax; i++){
        for(k = 0; k < n/circ_size; k++){
          if(array_element_in_range(col_row[circ_size*k+j], col_weight[circ_size*k+j], circ_size*i, circ_size*(i+1)-1)){
            double sum  = q0[circ_size*k+j];
            for(l = 0; l < col_weight[circ_size*k+j]; l++){
              sum += qin[circ_size*k+j][l];
            }
            //display_two_comp(sum - qin[circ_size*k+j][i]);
            //fprintf(Verilog_sim_f, "%08x (%d, %d, %d)\n", (int)((sum - qin[circ_size*k+j][var_neighbors_displayed[circ_size*k+j]++])*FRAC_LEVELS), j, i, k);
            fprintf(Verilog_sim_f, "%08x\n", (int)((sum - qin[circ_size*k+j][var_neighbors_displayed[circ_size*k+j]++])*FRAC_LEVELS));
          }
        }
      }
    }

    for (j = 0; j < m; j++) {

      int sgnprod = s[j];
      for (k = 0; k < row_weight[j]; k++){
        sgnprod ^= Sgntin[j][k];
      }

      for (k = 0; k < row_weight[j]; k++) {
        double min_msg = 200;
        for (l = 0; l < row_weight[j]; l++) {
          if(l != k){
            if(fabs(LogTanhtin[j][l]) < min_msg){
              min_msg = fabs(LogTanhtin[j][l]);
            }
          }
        }
        if(sgnprod != Sgntin[j][k]){
          *qin_row[j][k] = -min_msg;
        }
        else{
          *qin_row[j][k] = min_msg;
        }
      }
    }

    // For Verilog decoder debugging (for QC-code only)
    for(j = 0; j < circ_size; j++){
      for(k = 0; k < m/circ_size; k++){
        for(i = 0; i < row_weight[circ_size*k]; i++){
          //display_two_comp(*qin_row[circ_size*k+j][i]);
          //fprintf(Verilog_sim_f, "%08x (%d, %d, %d)\n", (int)(*qin_row[circ_size*k+j][i]*FRAC_LEVELS), j, k, i);
          fprintf(Verilog_sim_f, "%08x\n", (int)(*qin_row[circ_size*k+j][i]*FRAC_LEVELS));
        }
      }
    }
    
    for (i = 0; i < n; i++) {
      double sum = q0[i];
      for (j = 0; j < col_weight[i]; j++) {
        sum += qin[i][j];
      }
      sum = float_to_fix(sum);

      tmp_bit[i] = (sum < 0) ? 1 : 0;
    }
    
    enc(tmp_bit, tmp_s);
    i = HamDist(s, tmp_s, m);
    
    if (i == 0)           // nothing more can be done
      return 0;
  }

  return 1;
}

void initdec(char *s)
{
  int **row_N;
  int **col_N;
  int i, j, *count;
  FILE *fp = fopen(s, "rt");
  if (fp == NULL) {
    fprintf(stderr, "cannot open %s\n", s);
    exit(-2);
  }

  fscanf(fp, "%d%d", &n, &m);
  fscanf(fp, "%d%d", &cmax, &rmax);
  col_weight = malloc(sizeof(int) * n);
  for (i = 0; i < n; i++) {
    fscanf(fp, "%d", &col_weight[i]);
  }
  row_weight = malloc(sizeof(int) * m);
  for (j = 0; j < m; j++)
    fscanf(fp, "%d", &row_weight[j]);

  {//skip n lines
    for (i = 0; i < n; i++) {
      for (j = 0; j < cmax; j++)
        fscanf(fp, "%*d");
    }
  }

  count = malloc(sizeof(int) * n);
  memset(count, 0, sizeof(int) * n);
  qin = malloc2Ddouble(n, cmax);
  qin_row = malloc2Ddoublep(m, rmax);
  LogTanhtin     = malloc2Ddouble(m, rmax);
  LogTanhtin_row = malloc2Ddoublep(n, cmax);
  Sgntin     = malloc2Dint(m, rmax);
  Sgntin_row = malloc2Dintp(n, cmax);
  tmp_bit = malloc(sizeof(int) * n);
  tmp_s = malloc(sizeof(int) * m);

  var_neighbors_displayed = malloc(sizeof(int) * n);  // for the Verilog decoder debugging in dec()
  
  row_col = malloc2Dint(m, rmax);
  row_N   = malloc2Dint(m, rmax);
  col_row = malloc2Dint(n, cmax);
  col_N   = malloc2Dint(n, cmax);
  for (j = 0; j < m; j++) {
    for (i = 0; i < row_weight[j]; i++) {
      int v;
      fscanf(fp, "%d", &v);
      v--;
      row_col[j][i] = v;	// col address
      row_N[j][i] = count[v];	// vertical count of non-zero coef
      col_row[v][count[v]] = j;	// row address
      col_N[v][count[v]] = i;	// horizontal count of non-zero coef
      count[v]++;
      qin_row[j][i] = &qin[row_col[j][i]][row_N[j][i]];
    }
    // following block added on 02/05/2008 according to Mr. David Elkouss' comment
    for ( ; i < rmax; i++) {
      fscanf(fp, "%*d"); // skip the 0s (fillers)
    }
  }
  fclose(fp);

  for (i = 0; i < n; i++) {
    for (j = 0; j < col_weight[i]; j++) {
      LogTanhtin_row[i][j] = &LogTanhtin[col_row[i][j]][col_N[i][j]];
      Sgntin_row[i][j] =     &Sgntin    [col_row[i][j]][col_N[i][j]];
    }
  }

  free(count);
  free(*row_N);
  free( row_N);
  free(*col_N);
  free( col_N);
}

void test_code_min_sum_B_cw_noise_gen(char *code_filename, int iteration, int trials, double p_bsc, double *errors){
  INT_LEVELS = pow(2,INT-1);
  FRAC_LEVELS = pow(2,FRAC);

  srand(time(NULL));
  int i, j, dec_result;
  int *s, *x, *y;
  double *q0;
  
  initdec(code_filename);
  q0= malloc(sizeof(double) * n);
  s = malloc(sizeof(int) * m);  // syndrome
  x = malloc(sizeof(int) * n);  // source
  y = malloc(sizeof(int) * n);  // side information
  
  
  errors[0] = 0;
  errors[1] = 0;

  // Start Timer
  clock_t start = clock(), diff;

  for (i = 0; i < trials; i++) {
    for (j = 0; j < n; j++) {
      x[j] = rand() & 1;
      //x[j] = 0;
    }
    
    enc(x, s);
    bsc(x, y, p_bsc, q0);
    
    dec_result = dec(q0, s, iteration, x);
    
    if(dec_result){
      errors[0]++;
    } else {
      if(HamDist(tmp_bit, x, n) != 0) errors[1]++;
    }
  }

  // End Timer
  diff = clock() - start;
  int msec = diff * 1000 / CLOCKS_PER_SEC;
  //printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);
}

void test_code_min_sum_B(char *code_filename, int iteration, int *x, int *s, int *y, double *q0, double *errors, char *Verilog_sim_filename){
  INT_LEVELS = pow(2,INT-1);
  FRAC_LEVELS = pow(2,FRAC);

  srand(time(NULL));
  int i, dec_result;
  
  initdec(code_filename);
  
  errors[0] = 0;
  errors[1] = 0;

  Verilog_sim_f = fopen(Verilog_sim_filename, "wt");

  // Start Timer
  clock_t start = clock(), diff;
  
  dec_result = dec(q0, s, iteration, x);

  if(dec_result){
    errors[0]++;
  } else {
    if(HamDist(tmp_bit, x, n) != 0) errors[1]++;
  }

  fclose(Verilog_sim_f);

  // End Timer
  diff = clock() - start;
  int msec = diff * 1000 / CLOCKS_PER_SEC;
  //printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);
}


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    int i;

    char *code_filename;
    char *Verilog_sim_filename;
    int max_iter;
    int num_trials;
    double p_flip;
    int cw_noise_gen;

    int nvar;
    int nchk;
    double *x_in;
    double *s_in;
    double *y_in;
    double *q0;

    int *x;
    int *s;
    int *y;
    
    /* check for proper number of arguments */
    /* get the values of the inputs  */
    if(nrhs==4) {
        code_filename = mxArrayToString(prhs[0]);
        max_iter = mxGetScalar(prhs[1]);
        num_trials = mxGetScalar(prhs[2]);
        p_flip = mxGetScalar(prhs[3]);
        cw_noise_gen = 1;
    } else if(nrhs==9) {
        code_filename = mxArrayToString(prhs[0]);
        max_iter = mxGetScalar(prhs[1]);
        cw_noise_gen = 0;

        nvar = mxGetScalar(prhs[2]);
        nchk = mxGetScalar(prhs[3]);
        x_in = mxGetPr(prhs[4]);
        s_in = mxGetPr(prhs[5]);
        y_in = mxGetPr(prhs[6]);
        q0 = mxGetPr(prhs[7]);
        Verilog_sim_filename = mxArrayToString(prhs[8]);

        x = malloc(nvar*sizeof(int));
        s = malloc(nchk*sizeof(int));
        y = malloc(nvar*sizeof(int));

        for(i = 0; i < nvar; i++){
          x[i] = (int) x_in[i];
          y[i] = (int) y_in[i];
        }
        for(i = 0; i < nchk; i++){
          s[i] = (int) s_in[i];
        }
    } else{
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Number of inputs is wrong.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","Two output required.");
    }
        
    /* create the output matrix */
    plhs[0] = mxCreateNumericMatrix(1,2,mxDOUBLE_CLASS,mxREAL);
    
    /* get a pointer to the real data in the output matrix */
    double *outMatrix = mxGetPr(plhs[0]);

    /* call the computational routine */
    if(cw_noise_gen){
      test_code_min_sum_B_cw_noise_gen(code_filename, max_iter, num_trials, p_flip, outMatrix);
    } else{
      test_code_min_sum_B(code_filename, max_iter, x, s, y, q0, outMatrix, Verilog_sim_filename);
    }
}
