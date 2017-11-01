// GFq_LDPC_NTT.c, a GF(q>2) LDPC encoding/decoding simulator,
// powered by NTT (NumberTheoretic Transform, i.e., FFT over GF)
// Only allows GF(2^p), 2<=p<=8
// based on paper of M. C. Davey et al. "Low-Density...over GF(q)" June 1998
// (c) 2005-2006 by Seishi Takamura @ Stanford University / NTT (Nippon Telegraph and Telephone)
// Absolutely no warranty.
#include <mex.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <errno.h>
#include <time.h>


#define Q		13
#define Log2Q 4
#define half_Q_binary_ext  1<<(Log2Q-1)
#define GF_add(a, b) (a+b)%Q
#define GF_sub(a, b) (a-b)%Q


float **p_sent_given_rec_T;

int n, m;
int rmax, cmax;
int *row_weight, *col_weight;
int **row_col, **row_N;
int **col_row, **col_N;
int **H;
float ***logra;
float ***logqa;
float ***logracol;
float ***logqacol;
float **fQa;
int **isnegative;
int *tmp_z;
int *tmp_x;


int GF_mul(int a, int b)
{
  return (a*b)%Q;
}

// p_rec_given_sent[i][j] = P(i rec | j sent)
// p_sent_given_rec_T[i][j] = P(j sent | i rec)
void make_p_sent_given_rec_T(double** p_rec_given_sent){
  int i,j;
  float P_x = 1/(1.0*Q);
  float P_y;
  
  for(i = 0; i < Q; i++){
    P_y = 0;
    for(j = 0; j < Q; j++){
      P_y += (float) p_rec_given_sent[i][j]*P_x;
    }
    for(j = 0; j < Q; j++){
      p_sent_given_rec[i][j] = P_x*(float) p_rec_given_sent[i][j]/P_y;
    }
  }
}

// Modified from Takamura's NTT, by padding zeros to prime-sized
// data sequence, and taking NTT over the next binary power size
void ntt(float p[half_Q_binary_ext << 1])
{
  int b, factor = 1, rest;
  for (b = 0; b < Log2Q; b++) {
    for (rest = 0; rest < half_Q_binary_ext; rest++) {
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

void channel(int x[], int y[], double** p_rec_given_sent, float **logfna){
  int i, rand_select, rec_ind;
  float temp;

  for(i = 0; i < n; i++){
    temp = 0;
    rec_ind = 0;
    rand_select = rand()%101;
    while(temp <= rand_select){
      temp += (float) p_rec_given_sent[rec_ind++][x[i]];
    }
    assert(rec_ind > 0);
    y[i] = --rec_ind;
    logfna[i] = p_sent_given_rec_T[y[i]];
  }
}

int HamDist(int x[], int y[], int len)
{
  int i, j, sum = 0;
  for (i = 0; i < len; i++) {
    int xy = *x++ ^ *y++;
    for (j = 1; j < Q; j <<= 1)
      if ((xy & j) != 0) sum++;
  }
  return sum;
}

double mse(int x[], int y[], int len)
{
  int i;
  double sum = 0;
  for (i = 0; i < len; i++) {
    int s = *x++ - *y++;
    sum += (double)(s*s);
  }
  return sum / len;
}

// y[n] := x[n] + BSC_noise
// logfna[a][n]: log(Prob(x[n]=a | y[n]))
// can only be used for GF(2^b) case
void bsc(int x[], int y[], double p, int **logfna)
{
  int i, len = Log2Q * n;
  int modify = (int)(len * p + 0.5);
  int *err = malloc(sizeof(int) * n);
  double lp, l1p;
  p = modify / (double)len; // correct error probability
  //printf("m/n=%g, ", (double)m/n);
  //printf("BSC channel capacity(rate) = %g (bits)\n",
  //       -p*log2(p)-(1-p)*log2(1-p));
  memset(err, 0, sizeof(int) * n);
  // make sure p errors in err[]
  while (modify) {
    i = Rand() % len;
    if ((err[i/Log2Q] & (1<<(i%Log2Q))) != 0) continue;
    err[i/Log2Q] |= (1<<(i%Log2Q));
    modify--;
  }
  memcpy(y, x, sizeof(int) * n);
  for (i = 0; i < n; i++) {
    y[i] ^= err[i];
  }
  free(err);

  lp = log2(p);
  l1p = log2(1-p);
  for (i = 0; i < n; i++) {
    int a, j;
    for (a = 0; a < Q; a++) {
      double logprod = 0;
      for (j = 1; j < Q; j <<= 1) {
        if ((a&j) == (y[i]&j)) logprod += l1p;
        else logprod += lp;
      }
      logfna[i][a] = float2fix(logprod);
    }
  }
}

// y[n] := x[n] + Laplacian_noise
// logfna[n][a]: log(Prob(x[n]=a | y[n]))
// stddev: noise level
void lap(int x[], int y[], double stddev, int **logfna)
{
  int i;
  int count[Q];
  double _logfna[Q];
  double sum;
  memset(count, 0, sizeof(count));
  for (i = 0; i < n; i++) {
    double u2 = (Rand()+1) * (1.0 / (RandMax+1.0)); // uniform(0,1]
    int logu2 = (int)floor(log(u2) * stddev + .5);
    int a;
    sum = 0;
    if ((Rand() & 1) == 0) {
      y[i] = x[i] - logu2;
    } else {
      y[i] = x[i] + logu2;
    }
    // clipping
    if (y[i] < 0) y[i] = 0;
    else if (y[i] >= Q) y[i] = Q-1;
    //
    count[abs(x[i] - y[i])]++;
    //fna
    for (a = 0; a < Q; a++) {
      double dfna;
      //logfna[i][a] = -(abs(y[i] - a))/stddev;
      if (y[i] == a) _logfna[a] = -0.5/*0.2396*//stddev; //l((1-e(-0.5))*2)
      else _logfna[a] = -(abs(y[i] - a)-/*+*/0.0413)/stddev; //l(e(0.5)-e(-0.5))
      dfna = exp2(_logfna[a]);
      sum += dfna;
    }
    //normalize
    sum = log2(sum);
    for (a = 0; a < Q; a++) {
      logfna[i][a] = float2fix(_logfna[a] - sum);
    }
  }
  sum = 0;
  for (i = 0; i < Q; i++) {
    if (count[i]) sum += count[i] * log2(count[i]);
  }
  printf("m/n=%g, ", (double)m/n);
  printf("noise entropy = %g bits, rate = %g\n", -sum/n + log2(n),
         (-sum/n + log2(n)) / log2(Q));
  printf("PSNR = %g\n", 10*log((Q-1)*(Q-1)/mse(x,y,n)));
}

void enc(int x[], int s[])
{
  int i, j, xHmn;
  for (j = 0; j < m; j++) {
    int sum = 0;
    for (i = 0; i < row_weight[j]; i++) {
      xHmn = GF_mul(x[row_col[j][i]], H[j][i]);
      sum = GF_add(sum, xHmn);
    }
    s[j] = sum;
  }
}

// allocates float array[m][n]
float** malloc_2D_float(int m, int n){
  int i;
  float **pp = malloc(sizeof(float *) * m);
  float *p = malloc(sizeof(float) * m * n);
  if (pp == NULL || p == NULL) exit(-1);
  for (i = 0; i < m; i++) {
    pp[i] = p + n*i;
  }
  return pp;
}

// allocates float array[a][b][c]
float ***malloc_3D_float(int a, int b, int c)
{
  int i, j;
  float ***ppp= malloc(sizeof(float **) * a);
  float **pp  = malloc(sizeof(float *) * a * b);
  float *p    = malloc(sizeof(float) * a * b * c);
  if (ppp == NULL || pp == NULL || p == NULL) exit(-1);
  for (j = 0; j < a; j++) {
    for (i = 0; i < b; i++) {
      pp[i+b*j] = p + c*i + c*b*j;
    }
    ppp[j] = pp + b*j;
  }
  return ppp;
}

// allocates array[a][b] (float pointer)
float ***malloc_2D_float_p(int a, int b)
{
  int j;
  float ***ppp= malloc(sizeof(float **) * a);
  float **pp  = malloc(sizeof(float *) * a * b);
  if (ppp == NULL || pp == NULL) exit(-1);
  for (j = 0; j < a; j++) {
    ppp[j] = pp + b*j;
  }
  return ppp;
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

// Sum Product Decoder
// logfna: prior
// z: syndrome
// loop_max: max iteration
// x[]: original signal (just for reference)
int dec(int **logfna, int z[], int loop_max)
{
  int a, i, j, k, loop;
  float normalizing_val;

  for (i = 0; i < n; i++) {
    for (k = 0; k < col_weight[i]; k++) {
      memcpy(logqacol[i][k], logfna[i], sizeof(float)*Q);
    }
  }

  for (loop = 0; loop < loop_max; loop++) {
    for (j = 0; j < m; j++) {
      int row_weightj = row_weight[j];
      float logprodqa[half_Q_binary_ext << 1];
      //int sgnsum[Q];
      //memset(sgnsum, 0, sizeof(sgnsum));
      memset(logprodqa, 1.0, sizeof(logprodqa));
      //memset(*isnegative, 0.0, rmax * Q * sizeof(int));
      for (k = 0; k < row_weightj; k++) {
        for (a = 0; a < Q; a++) {
          fQa[k][GF_mul(a, H[j][k])] = logqa[j][k][a];
        }
        memset(fQa[k]+Q, 0, (2*half_Q_binary_ext-Q)*sizeof(float));    // Pad Zeros
        
        ntt(fQa[k]);

        for (a = 0; a < half_Q_binary_ext << 1; a++) {
          logprodqa[a] *= fQa[k][a];
        }
      }
      for (k = 0; k < row_weightj; k++) {
        float fRa[half_Q_binary_ext << 1];
        for (a = 0; a < half_Q_binary_ext << 1; a++) {
          fRa[a] = logprodqa[a]/fQa[k][a];
        }
        ntt(fRa);
        normalizing_val = 0;
        for (a = 0; a < Q; a++) {
          logra[j][k][a] = fRa[GF_sub(z[j], GF_mul(a, H[j][k]))];
          normalizing_val += logra[j][k][a];
        }
        for (a = 0; a < Q; a++) {
          logra[j][k][a] = logra[j][k][a]/normalizing_val;
        }
      }
    }
    // update qa & tentative decode
    for (i = 0; i < n; i++) {
      int col = col_weight[i];
      float argmaxa = -1;
      float maxp = -(1<<30);
      for (a = 0; a < Q; a++) {
        int logprod = logfna[i][a];
        for (k = 0; k < col; k++) {
          logprod += logracol[i][k][a];
        }
        for (k = 0; k < col; k++) {
          logqacol[i][k][a] = logprod - logracol[i][k][a]; //taka product except k
        }
        if (logprod > maxp) {
          maxp = logprod;
          argmaxa = a;
        }
      }
      tmp_x[i] = argmaxa;
    }
    //normalize qa
    for (j = 0; j < m; j++) {
      for (k = 0; k < row_weight[j]; k++) {
        float *qajk = logqa[j][k];
        float sum = 0;
        for (a = 0; a < Q; a++) {
          sum += *qajk++;
        }/*
        if (sum <= 0) {//all zero happenes
          //printf("%d,%d:sum%d ", j, k, sum);
          sum = Q;//taka051207
        }*/
        qajk = logqa[j][k];
        for (a = 0; a < Q; a++) {
          *qajk = *qajk/sum;
          qajk++;
        }
      }
    }

    enc(tmp_x, tmp_z);
    {
      int dist = HamDist(z, tmp_z, m);
      //printf("HamDist(s,synd(x^))=%d\n", dist);
      if (dist == 0)           // nothing more can be done
        return 0;
    }
  }

  return 1;
}

//http://www.inference.phy.cam.ac.uk/mackay/codes/alist.html
void initdec(char *s)
{
  int i, j, q, *count;
  char buf[BUFSIZ];
  FILE *fp = fopen(s, "rt");
  if (fp == NULL) {
    fprintf(stderr, "cannot open %s\n", s);
    exit(-2);
  }
  fgets(buf, sizeof(buf), fp);
  if (sscanf(buf, "%d%d%d", &n, &m, &q) == 2) {
    fprintf(stderr, "Warning! A binary matrix seems to be specified.\n");
    fprintf(stderr, "I'll use it anyways.\n");
  }
  if (q > Q) exit(-1);
  if (q != Q) {
    fprintf(stderr, "Warning! GF order(%d) does not match with matrix file(%d)\n", Q, q);
    fprintf(stderr, "I'll use it anyways.\n");
  }
  fscanf(fp, "%d%d", &cmax, &rmax);
  col_weight = malloc(sizeof(int) * n);
  for (i = 0; i < n; i++) {
    fscanf(fp, "%d", &col_weight[i]);
  }
  row_weight = malloc(sizeof(int) * m);
  for (j = 0; j < m; j++)
    fscanf(fp, "%d", &row_weight[j]);

  count = malloc(sizeof(int) * n);
  memset(count, 0, sizeof(int) * n);

  row_col = malloc2Dint(m, rmax);
  row_N   = malloc2Dint(m, rmax);
  col_row = malloc2Dint(n, cmax);
  col_N   = malloc2Dint(n, cmax);
  H            = malloc2Dint(m, rmax); // Hmn
  fQa          = malloc_2D_float(rmax, half_Q_binary_ext << 1);
  isnegative   = malloc2Dint(rmax, Q);

  {//skip n lines
    for (i = 0; i < n; i++) {
      for (j = 0; j < cmax; j++)
        fscanf(fp, "%*d%*d");
    }
  }

  for (j = 0; j < m; j++) {
    for (i = 0; i < row_weight[j]; i++) {
      int v;
      fscanf(fp, "%d%d", &v, &H[j][i]);
      v--;
      // v (at this point) is the absolute index of the column in the H matrix

      row_col[j][i] = v;  // row_col[j][i] is the absolute index of the i-th variable neighor of the j-th check
      row_N[j][i] = count[v];   // row_N[j][i] the relative index of the j-th check neighbor to the i-th variable
      col_row[v][count[v]] = j;   // col_row[j][i] the absolute index of the i-th check neighbor to the j-th variable
      col_N[v][count[v]] = i;   // col_N[j][i] relative index of the j-th variable in its i-th check
      count[v]++;
    }

    for (; i < rmax; i++) {//skip fillers
      int a,b;
      fscanf(fp, "%d%d", &a, &b);
      if (a!=0 || b!=0) {
        printf("error at row %d, %d %d\n", i, a, b);
        exit(-1);
      }
    }
  }
  free(count);
  fclose(fp);
  tmp_z = malloc(sizeof(int) * m);
  tmp_x = malloc(sizeof(int) * n);

  logra  = malloc_3D_float(m, rmax, Q);
  logqa  = malloc_3D_float(m, rmax, Q);
  logracol=malloc_2D_float_p(n, cmax);
  logqacol=malloc_2D_float_p(n, cmax);
  for (i = 0; i < n; i++) {
    int k;
    int *cri = col_row[i];  // pointer to neighbor check indices
    int *cNi = col_N[i];  // pointer to relative indices to the variable in its checks
    int col = col_weight[i];
    for (k = 0; k < col; k++) {
      logracol[i][k] = logra[cri[k]][cNi[k]];
      logqacol[i][k] = logqa[cri[k]][cNi[k]];
    }
  }
}

void test_NB_LDPC_p(int iteration, int num_trials, double** p_rec_given_sent, double *errors)
{
  int i, j, *s, *x, *y, dec_result;
  float **logfna;


  initdec("bazar_GF_13.txt");
  p_sent_given_rec_T = malloc_2D_float(Q, Q);
  s = malloc(sizeof(int) * m);  // syndrome
  x = malloc(sizeof(int) * n);  // source
  y = malloc(sizeof(int) * n);  // side information
  logfna = malloc_2D_float(n, Q);
  make_p_sent_given_rec_T(p_rec_given_sent);

  errors[0] = 0;
  errors[1] = 0;

  // Start Timer
  clock_t start = clock(), diff;

  // iterate experiments
  for (j = 1; j <= num_trials; j++) {
    srand(j);
    // Generate random data
    for (i = 0; i < n; i++)
      x[i] = rand() % Q;

    enc(x, s);
    channel(x, y, p_rec_given_sent, logfna);
    dec_result = dec(logfna, s, iteration);
    
    if(dec_result){
      errors[0]++;
    } else {
      if(HamDist(tmp_x, x, n) != 0) errors[1]++;
    }
  }

  // End Timer
  diff = clock() - start;
  int msec = diff * 1000 / CLOCKS_PER_SEC;
  printf("Time taken %d seconds %d milliseconds", msec/1000, msec%1000);
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double **conf_mat;              /* input scalar */
    int num_trials;              /* input scalar */
    int max_iter;              /* input scalar */
    
    /* check for proper number of arguments */
    if(nrhs!=3) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Three inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
    }
    
    /* get the values of the inputs  */
    max_iter = mxGetScalar(prhs[0]);
    num_trials = mxGetScalar(prhs[1]);
    conf_mat = mxGetPr(prhs[2]);
    
    /* create the output matrix */
    plhs[0] = mxCreateNumericMatrix(1,2,mxDOUBLE_CLASS,mxREAL);
    
    /* get a pointer to the real data in the output matrix */
    double *outMatrix = mxGetPr(plhs[0]);

    /* call the computational routine */
    test_NB_LDPC_p(max_iter, num_trials, conf_mat, outMatrix);
}
