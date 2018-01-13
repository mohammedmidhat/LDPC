#include "mex.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>


#define Q 16
#define page_size 18336
#define QLC   4
#define CW_per_page   2
#define num_sym_per_bit		8


int grey_code_inv[16] = {6,5,7,14,9,12,8,13,3,4,2,15,10,11,1,0};
int sym_top_1[num_sym_per_bit] = {0,1,2,3,4,10,11,15};
int sym_up_1[num_sym_per_bit] = {0,1,8,9,10,11,12,13};
int sym_mid_1[num_sym_per_bit] = {0,1,2,7,8,13,14,15};
int sym_low_1[num_sym_per_bit] = {0,4,5,11,12,13,14,15};

int n, m;
int rmax, cmax;
int *row_weight, *col_weight;
int **row_col;
double **p_sent_given_rec_T;

// Lower page data stored first
char data_read[8*page_size];
char data_written[8*page_size];

double atanh2(double x)
{
  return log((1.0 + x) / (1.0 - x));  // returns 2*atanh(x)
}
double logtanh2(double x)
{
  return log(tanh(fabs(x*0.5)));  // returns log tanh |x|
}

#define INT	8/*8*/              // int part
#define DECI	14/*13*/              // fraction part
#define FMUL	(1<<DECI)       // multiplier
#define PREC	(1.0/FMUL)      // precision
#define LEVELS	(1<<(INT+DECI))
static int flogtanh[LEVELS];
static int fgallag[LEVELS];

int float2fix(double x)
{
  if (x >= 0) {
    return (int)(x * FMUL + 0.5);
  } else {
    return -(int)((-x) * FMUL + 0.5);
  }
}

unsigned int float2fixu(double x)
{
  return (unsigned int)(x * FMUL + 0.5);
}

#define fix2float(x) ((x)*PREC)

void inittab(void)
{
  int i = 1;
  double right = logtanh2(fix2float(i) - 0.5*PREC);
  flogtanh[0] = -FMUL*14;
  for ( ; i < LEVELS; i++) {
    double d = fix2float(i);
    double left = logtanh2(d+0.5*PREC);
    flogtanh[i] = float2fix((4*logtanh2(d)+right+left) / 6.0);
    right = left;
  }

  i = 1;
  fgallag[0] = FMUL*14;
  right = atanh2(exp(fix2float(-i) - 0.5*PREC));
  for ( ; i < LEVELS; i++) {
    double d = fix2float(-i);
    double expd = atanh2(exp(d));
    double left = atanh2(exp(d+0.5*PREC));
    fgallag[i] = float2fix((4*expd+right+left) / 6.0);
    right = left;
  }
}

int Flogtanh(int x)
{
  assert(x>=0);//if (x < 0) return 0;
  if (x >= LEVELS)
    return 0;
  return flogtanh[x];
}

int Fgallag(int x)
{
  assert(x <= 0);//  if (x >= 0) return -FMUL*14; //-115000
  if (x <= -LEVELS)
    return 0;
  return fgallag[-x];
}

int HamDist(int *x, int *y, int len)
{
  int i, sum = 0;
  for (i = 0; i < len; i++) {
    if (*x++ != *y++) sum++;
  }
  return sum;
}

int bsc(int x[], int y[], double p, int q0[])
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
    q0[i] = float2fix(d);
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

int **qin, ***qin_row;
int **LogTanhtin, ***LogTanhtin_row;
int **Sgntin, ***Sgntin_row;
int *tmp_bit;
int *tmp_s;
int loop;

int dec(int q0[], int s[], int loop_max)
{
  int i, j, k;
  int iir, prev = 999999, nodecr = 0;

  memset(*qin, 0, n * cmax * sizeof(int));

  for (loop = 0; loop < loop_max; loop++) {
    for (i = 0; i < n; i++) {
      int sum = q0[i];
      for (j = 0; j < col_weight[i]; j++)
        sum += qin[i][j];
      for (j = 0; j < col_weight[i]; j++) {
        int qout = sum - qin[i][j];
        if (qout < 0) {
          *LogTanhtin_row[i][j] = Flogtanh(-qout);
          *Sgntin_row[i][j] = 1;
        } else {
          *LogTanhtin_row[i][j] = Flogtanh(qout);
          *Sgntin_row[i][j] = 0;
        }
        //printf("v_msg_%d_%d = %d",i,j,qout);
      }
    }

    for (j = 0; j < m; j++) {
      int sgnprod = s[j];
      int logprod = 0;
      for (k = 0; k < row_weight[j]; k++) {
        logprod += LogTanhtin[j][k];
        sgnprod ^= Sgntin[j][k];
      }

      for (k = 0; k < row_weight[j]; k++) {
        int tout = Fgallag(logprod - LogTanhtin[j][k]);
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
      int sum = q0[i];
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
  qin = malloc2Dint(n, cmax);
  qin_row = malloc2Dintp(m, rmax);
  LogTanhtin     = malloc2Dint(m, rmax);
  LogTanhtin_row = malloc2Dintp(n, cmax);
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
    /*for ( ; i < rmax; i++) {
      fscanf(fp, "%*d"); // skip the 0s (fillers)
    }*/
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

// p_rec_given_sent[i][j] = P(i rec | j sent)
// p_sent_given_rec_T[i][j] = P(j sent | i rec)
void make_p_sent_given_rec_T(double* p_rec_given_sent, int num_reads){
  int i,j;
  double P_x = 1/(1.0*Q);
  double P_y;

  for(i = 0; i < Q*num_reads; i++){
    P_y = 0;
    for(j = 0; j < Q; j++){
      P_y += p_rec_given_sent[i*Q+j]*P_x;
    }
    if(P_y){
      for(j = 0; j < Q; j++){
        p_sent_given_rec_T[i][j] = P_x*p_rec_given_sent[i*Q+j]/P_y;
      }
    } else{
      for(j = 0; j < Q; j++){
        p_sent_given_rec_T[i][j] = 0;
      }
    }
  }
}

void get_bits_in_symbol(char symbol, int x[], int symbol_ind){
  int lp_bit, mp_bit, up_bit, tp_bit;

  lp_bit = symbol & 1;
  mp_bit = (symbol >> 1) & 1;
  up_bit = (symbol >> 2) & 1;
  tp_bit = (symbol >> 3) & 1;

  x[4*symbol_ind] = lp_bit;
  x[4*symbol_ind+1] = mp_bit;
  x[4*symbol_ind+2] = up_bit;
  x[4*symbol_ind+3] = tp_bit;
}

void channel(int x[], int y[], double* p_rec_given_sent, int q0[], int num_reads){
  int i, symbol, rand_select, rec_ind;
  double temp;

  for(i = 0; i < n/4; i++){
    symbol = grey_code_inv[(x[4*i+3] << 3) + (x[4*i+2] << 2) + (x[4*i+1] << 1) + (x[4*i])];

    temp = 0;
    rec_ind = 0;
    rand_select = rand()%101;

    if(!rand_select){
      while(!p_rec_given_sent[(rec_ind++)*Q+x[i]]);
    }

    while(temp < rand_select && rec_ind != Q*num_reads-2){
      temp += 100*p_rec_given_sent[(rec_ind++)*Q+x[i]];
    }

    symbol = --rec_ind;

    assign_llr_one_sym(symbol, i, q0);
  }
}

void assign_llr_one_sym(int sym, int sym_ind, int q0[]){
  double Pr_1, llr;
  int j;

  // On lower page
    Pr_1 = 0;
    for(j = 0; j < num_sym_per_bit; j++){
      Pr_1 += p_sent_given_rec_T[sym][sym_low_1[j]];
    }
    if(Pr_1 == 1){
      llr = -100;
    } else if(Pr_1 == 0){
      llr = 100;
    } else{
      llr = log((1.0 - Pr_1) / Pr_1);
    }
    q0[4*sym_ind] = float2fix(llr);
    // On middle page
    Pr_1 = 0;
    for(j = 0; j < num_sym_per_bit; j++){
      Pr_1 += p_sent_given_rec_T[sym][sym_mid_1[j]];
    }
    if(Pr_1 == 1){
      llr = -100;
    } else if(Pr_1 == 0){
      llr = 100;
    } else{
      llr = log((1.0 - Pr_1) / Pr_1);
    }
    q0[4*sym_ind+1] = float2fix(llr);
    // On upper page
    Pr_1 = 0;
    for(j = 0; j < num_sym_per_bit; j++){
      Pr_1 += p_sent_given_rec_T[sym][sym_up_1[j]];
    }
    if(Pr_1 == 1){
      llr = -100;
    } else if(Pr_1 == 0){
      llr = 100;
    } else{
      llr = log((1.0 - Pr_1) / Pr_1);
    }
    q0[4*sym_ind+2] = float2fix(llr);
    // On top page
    Pr_1 = 0;
    for(j = 0; j < num_sym_per_bit; j++){
      Pr_1 += p_sent_given_rec_T[sym][sym_top_1[j]];
    }
    if(Pr_1 == 1){
      llr = -100;
    } else if(Pr_1 == 0){
      llr = 100;
    } else{
      llr = log((1.0 - Pr_1) / Pr_1);
    }
    q0[4*sym_ind+3] = float2fix(llr);
}

// Works only if n is multiple of 4 (a single QLC cell is not shared
// between 2 CWs)
void assign_llr(int y[], int q0[]){
  int i, j;
  char symbol;
  double Pr_1, llr;

  for(i = 0; i < n/4; i++){
    symbol = y[i];

    assign_llr_one_sym(symbol, i, q0);
  }
}

void test_code_B_MSDP(int iteration, int num_trials, int num_reads, int decode_mode, double* p_rec_given_sent, double *errors){
  int i, j, dec_result, *iterations, *s, *x, *y, *q0;
  int CW_per_page_fetched;
  char symbol;
  
  inittab();
  
  initdec("16383.2130.3.txt");
  p_sent_given_rec_T = malloc2Ddouble(Q*num_reads, Q);
  q0= malloc(sizeof(int) * n);
  s = malloc(sizeof(int) * m);  // syndrome
  x = malloc(sizeof(int) * n);  // source
  y = malloc(sizeof(int) * n);  // side information
  
  make_p_sent_given_rec_T(p_rec_given_sent, num_reads);

  errors[0] = 0;
  errors[1] = 0;

  // Start Timer
  clock_t start = clock(), diff;

  if(decode_mode){
    FILE *read_data = fopen("snowbird_sym.bin", "rb");
    FILE *written_data = fopen("snowbird_sym.bin", "rb");

    while(num_trials){
      fread(data_read, 1, 8*page_size, read_data);
      fread(data_written, 1, 8*page_size, written_data);
      CW_per_page_fetched = 0;

      while(CW_per_page_fetched != CW_per_page){
        for(i = 0; i < n/4; i++){
            symbol = data_written[CW_per_page_fetched*n/4+i];
            get_bits_in_symbol(symbol, x, i);
            symbol = data_read[CW_per_page_fetched*n/4+i];
            y[i] = symbol;
        }

        enc(x, s);
        assign_llr(y, q0);
    
        dec_result = dec(q0, s, iteration);
    
        if(dec_result){
            errors[0]++;
        } else {
            if(HamDist(tmp_bit, x, n) != 0) errors[1]++;
        }
        CW_per_page_fetched++;
        num_trials--;
      }
    }

    fclose(read_data);
    fclose(written_data);
  } else{
    for (j = 1; j <= num_trials; j++){
      srand(j);
      // Generate random data
      for (i = 0; i < n; i++){
        x[i] = rand() % 2;
      }

      enc(x, s);
      channel(x, y, p_rec_given_sent, q0, num_reads);

      dec_result = dec(q0, s, iteration);

      if(dec_result){
          errors[0]++;
      } else {
          if(HamDist(tmp_bit, x, n) != 0) errors[1]++;
      }
    }
  }
  
  // End Timer
  diff = clock() - start;
  int msec = diff * 1000 / CLOCKS_PER_SEC;
  printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);

  /*printf("Undetected Errors = %f\n", errors[1]);
  printf("Errors = %f\n", errors[0]);*/
}


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *conf_mat;              // input scalar
    int num_trials;              // input scalar
    int max_iter;              // input scalar
    int num_reads;               // input scalar
    int decode_mode;          // input scalar
    
    /* check for proper number of arguments */
    if(nrhs!=5) {
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
    conf_mat = mxGetPr(prhs[4]);
    
    /* create the output matrix */
    plhs[0] = mxCreateNumericMatrix(1,2,mxDOUBLE_CLASS,mxREAL);
    
    /* get a pointer to the real data in the output matrix */
    double *outMatrix = mxGetPr(plhs[0]);

    /* call the computational routine */
    test_code_B_MSDP(max_iter, num_trials, num_reads, decode_mode, conf_mat, outMatrix);
}
