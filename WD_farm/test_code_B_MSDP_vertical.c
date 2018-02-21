#include "mex.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>


#define Q 16


int n, m;
int rmax, cmax;
int *row_weight, *col_weight;
int **row_col;


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
    q0[i] = d;
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

int dec(double q0[], int s[], int loop_max)
{
  int i, j, k, loop;

  memset(*qin, 0, n * cmax * sizeof(double));

  for (loop = 0; loop < loop_max; loop++) {
    for (i = 0; i < n; i++) {
      double sum = q0[i];
      for (j = 0; j < col_weight[i]; j++)
        sum += qin[i][j];
      for (j = 0; j < col_weight[i]; j++) {
        double qout = sum - qin[i][j];
        if (qout < 0) {
          *LogTanhtin_row[i][j] = logtanh2(-qout);
          *Sgntin_row[i][j] = 1;
        } else {
          *LogTanhtin_row[i][j] = logtanh2(qout);
          *Sgntin_row[i][j] = 0;
        }
        //printf("v_msg_%d_%d = %d",i,j,qout);
      }
    }

    for (j = 0; j < m; j++) {
      int sgnprod = s[j];
      double logprod = 0;
      for (k = 0; k < row_weight[j]; k++) {
        logprod += LogTanhtin[j][k];
        sgnprod ^= Sgntin[j][k];
      }

      for (k = 0; k < row_weight[j]; k++) {
        double tout = atanh2(exp(logprod - LogTanhtin[j][k]));
        if(sgnprod != Sgntin[j][k]){
          *qin_row[j][k] = -tout;
          //printf("c_msg_%d_%d = %d",j,k,-tout);
        } else{
          *qin_row[j][k] = tout;
          //printf("c_msg_%d_%d = %d",j,k,tout);
        }
      }
    }

    for (i = 0; i < n; i++) {
      double sum = q0[i];
      for (j = 0; j < col_weight[i]; j++) {
        sum += qin[i][j];
      }

      tmp_bit[i] = (sum < 0) ? 1 : 0;
    }
    //printf("%2d:HamDist(x)=%d\n ", loop+1, HamDist(x, tmp_bit, n));

    enc(tmp_bit, tmp_s);
    i = HamDist(s, tmp_s, m);
    //printf("HamDist(s,synd(x^))=%d\n", i);
    if (i == 0)           // nothing more can be done
      return 0;
  }
  return 1;
}

void initdec(char *s)
{
  int **row_N;
  int **col_row, **col_N;
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
      for (j = 0; j < col_weight[i]; j++)
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

  row_col = malloc2Dint(m, rmax);
  row_N   = malloc2Dint(m, rmax);
  col_row = malloc2Dint(n, cmax);
  col_N   = malloc2Dint(n, cmax);
  for (j = 0; j < m; j++) {
    for (i = 0; i < row_weight[j]; i++) {
      int v;
      fscanf(fp, "%d", &v);
      v--;
      //if(j == 0){
      //  printf("%d\n", v);
      //}
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
  free(*col_row);
  free( col_row);
  free(*col_N);
  free( col_N);
}

void channel(int x[], int y[], double q0[], int num_reads, double lp_p, double mp_p, double up_p, double tp_p, FILE *debug_file){
  int i, rand_select;
  
  for(i = 0; i < n/4; i++){    
    rand_select = rand()%100;
    if(lp_p*100 > rand_select){
      y[4*i] = x[4*i]^1;
    }
    q0[4*i] = (1 - 2 * y[4*i]) * log((1.0 - lp_p) / lp_p);

    rand_select = rand()%100;
    if(mp_p*100 > rand_select){
      y[4*i+1] = x[4*i+1]^1;
    }
    q0[4*i+1] = (1 - 2 * y[4*i+1]) * log((1.0 - lp_p) / lp_p);

    rand_select = rand()%100;
    if(up_p*100 > rand_select){
      y[4*i+2] = x[4*i+2]^1;
    }
    q0[4*i+2] = (1 - 2 * y[4*i+2]) * log((1.0 - lp_p) / lp_p);

    rand_select = rand()%100;
    if(tp_p*100 > rand_select){
      y[4*i+3] = x[4*i+3]^1;
    }
    q0[4*i+3] = (1 - 2 * y[4*i+3]) * log((1.0 - lp_p) / lp_p);
  }
}

void test_code_B_MSDP_vertical(int iteration, int num_trials, int num_reads, int decode_mode, double lp_p, double mp_p, double up_p, double tp_p, double *errors){
  srand(time(NULL));
  int i, j, k, dec_result, *iterations, *s, *x, *y;
  double *q0;

  
  initdec("peg_16000_3_0.9.txt");
  q0= malloc(sizeof(double) * n);
  s = malloc(sizeof(int) * m);  // syndrome
  x = malloc(sizeof(int) * n);  // source
  y = malloc(sizeof(int) * n);  // side information
  

  errors[0] = 0;
  errors[1] = 0;

  // Start Timer
  clock_t start = clock(), diff;

  FILE *debug_f = fopen("debug.txt", "w");
  for (j = 1; j <= num_trials; j++){
    // Generate random data
    for (i = 0; i < n; i++){
      x[i] = rand() % 2;
    }

    enc(x, s);
    /*for(i = 0; i < n; i++){
      fprintf(debug_f, "%d ", x[i]);
    }
    fprintf(debug_f, "\n");*/
    channel(x, y, q0, num_reads, lp_p, mp_p, up_p, tp_p, debug_f);
    /*for(i = 0; i < n; i++){
      fprintf(debug_f, "%.3f ", q0[i]);
    }
    fprintf(debug_f, "\n");*/
    
    dec_result = dec(q0, s, iteration);

    if(dec_result){
        errors[0]++;
    } else {
        if(HamDist(tmp_bit, x, n) != 0) errors[1]++;
    }
  }
  fclose(debug_f);
  
  // End Timer
  diff = clock() - start;
  int msec = diff * 1000 / CLOCKS_PER_SEC;
  //printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);

  /*printf("Undetected Errors = %f\n", errors[1]);
  printf("Errors = %f\n", errors[0]);*/
}


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double lp_p;              // input scalar
    double mp_p;              // input scalar
    double up_p;              // input scalar
    double tp_p;              // input scalar
    int num_trials;              // input scalar
    int max_iter;              // input scalar
    int num_reads;               // input scalar
    int decode_mode;          // input scalar
    
    /* check for proper number of arguments */
    if(nrhs!=8) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Five inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","Two output required.");
    }
    
    /* get the values of the inputs  */
    max_iter = mxGetScalar(prhs[0]);
    num_trials = mxGetScalar(prhs[1]);
    num_reads = mxGetScalar(prhs[2]);
    decode_mode = mxGetScalar(prhs[3]);
    lp_p = mxGetScalar(prhs[4]);
    mp_p = mxGetScalar(prhs[5]);
    up_p = mxGetScalar(prhs[6]);
    tp_p = mxGetScalar(prhs[7]);
    
    /* create the output matrix */
    plhs[0] = mxCreateNumericMatrix(1,2,mxDOUBLE_CLASS,mxREAL);
    
    /* get a pointer to the real data in the output matrix */
    double *outMatrix = mxGetPr(plhs[0]);

    /* call the computational routine */
    test_code_B_MSDP_vertical(max_iter, num_trials, num_reads, decode_mode, lp_p, mp_p, up_p, tp_p, outMatrix);
}
