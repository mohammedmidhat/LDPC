// Binary LDPC decoder, called in from a MATLAB script


#include "mex.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>

int n, m;
int rmax, cmax;
int *row_weight, *col_weight;
int **row_col;

//int err_count, undet_err_count;

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

int **qin, ***qin_row;
int **LogTanhtin, ***LogTanhtin_row;
int **Sgntin, ***Sgntin_row;
int *tmp_bit;
int *tmp_s;


int dec(int q0[], int s[], int loop_max, int x[])
{
  int i, j, k, loop;
  
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
  //printf("1\n");
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
      //printf("2\n");
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

void test_code_B(int iteration, int trials, double p_bsc, double *errors){
  srand(time(NULL));
  int i, j, dec_result, *iterations, *s, *x, *y, *q0;
  int CW_per_page_fetched;
  
  inittab();
  
  initdec("peg_16000_3_0.9.txt");
  q0= malloc(sizeof(int) * n);
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
  printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);

  /*printf("Undetected Errors = %f\n", errors[1]);
  printf("Errors = %f\n", errors[0]);*/
}


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    int max_iter;
    int num_trials;
    double p_flip;
    
    /* check for proper number of arguments */
    if(nrhs!=3) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Three inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","Two output required.");
    }
    
    /* get the values of the inputs  */
    max_iter = mxGetScalar(prhs[0]);
    num_trials = mxGetScalar(prhs[1]);
    p_flip = mxGetScalar(prhs[2]);
    
    /* create the output matrix */
    plhs[0] = mxCreateNumericMatrix(1,2,mxDOUBLE_CLASS,mxREAL);
    
    /* get a pointer to the real data in the output matrix */
    double *outMatrix = mxGetPr(plhs[0]);

    /* call the computational routine */
    test_code_B(max_iter, num_trials, p_flip, outMatrix);
}
