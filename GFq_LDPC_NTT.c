// GFq_LDPC_NTT.c, a GF(q>2) LDPC encoding/decoding simulator,
// powered by NTT (NumberTheoretic Transform, i.e., FFT over GF)
// Only allows GF(2^p), 2<=p<=8
// based on paper of M. C. Davey et al. "Low-Density...over GF(q)" June 1998
// (c) 2005-2006 by Seishi Takamura @ Stanford University / NTT (Nippon Telegraph and Telephone)
// Absolutely no warranty.
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <errno.h>
#include <time.h>

#ifndef Log2Q
  #define Log2Q 5               // GF(2^5)
#endif
#if Log2Q < 1 || Log2Q > 8 
  #error "Log2Q must be 1..8"
#endif

#define Q (1<<Log2Q)      // GF(Q)

#if Q==4
const int logq[4] = {0,0,1,2};
const int expq[3] = {1,2,3};
#elif Q==8
const int logq[8] = {0,0,1,3,2,6,4,5};
const int expq[7] = {1,2,4,3,6,7,5};
#elif Q==16
const int logq[16] = {0,0,1,4,2,8,5,10,3,14,9,7,6,13,11,12};
const int expq[15] = {1,2,4,8,3,6,12,11,5,10,7,14,15,13,9};
#elif Q==32
const int logq[32] = {0,0,1,18,2,5,19,11,3,29,6,27,20,8,12,23,4,
  10,30,17,7,22,28,26,21,25,9,16,13,14,24,15};
const int expq[31] = {1,2,4,8,16,5,10,20,13,26,17,7,14,28,29,31,
  27,19,3,6,12,24,21,15,30,25,23,11,22,9,18};
#elif Q==64
const int logq[64] = {0,0,1,6,2,12,7,26,3,32,13,35,8,48,27,18,4,24,
  33,16,14,52,36,54,9,45,49,38,28,41,19,56,5,62,
  25,11,34,31,17,47,15,23,53,51,37,44,55,40,10,
  61,46,30,50,22,39,43,29,60,42,21,20,59,57,58};
const int expq[63] = {1,2,4,8,16,32,3,6,12,24,48,35,5,10,20,40,19,
  38,15,30,60,59,53,41,17,34,7,14,28,56,51,37,
  9,18,36,11,22,44,27,54,47,29,58,55,45,25,50,
  39,13,26,52,43,21,42,23,46,31,62,63,61,57,49,33};
#elif Q==128
const int logq[128] = {0,0,1,31,2,62,32,103,3,7,63,15,33,84,104,
  93, 4,124,8,121,64,79,16,115,34,11,85,38,105,46,94,51,
  5,82,125,60,9,44,122,77,65,67,80,42,17,69,116,23,35,118,
  12,28,86,25,39,57,106,19,47,89,95,71,52,110,6,14,83,92,126,
  30,61,102,10,37,45,50,123,120,78,114,66,41,68,22,81,59,43,76,
  18,88,70,109,117,27,24,56,36,49,119,113,13,91,29,101,87,108,
  26,55,40,21,58,75,107,54,20,74,48,112,90,100,96,97,72,98,53,73,111,99};
const int expq[127] = {1,2,4,8,16,32,64,9,18,36,72,25,50,100,65,11,
  22,44,88,57,114,109,83,47,94,53,106,93,51,102,69,3,6,12,24,
  48,96,73,27,54,108,81,43,86,37,74,29,58,116,97,75,31,62,124,
  113,107,95,55,110,85,35,70,5,10,20,40,80,41,82,45,90,61,122,
  125,115,111,87,39,78,21,42,84,33,66,13,26,52,104,89,59,118,101,
  67,15,30,60,120,121,123,127,119,103,71,7,14,28,56,112,105,91,63,
  126,117,99,79,23,46,92,49,98,77,19,38,76,17,34,68};
#elif Q==256
const int logq[256] = {0,0,1,25,2,50,26,198,3,223,51,238,27,104,199,75,4,100,
  224,14,52,141,239,129,28,193,105,248,200,8,76,113,5,138,101,47,225,
  36,15,33,53,147,142,218,240,18,130,69,29,181,194,125,106,39,249,185,
  201,154,9,120,77,228,114,166,6,191,139,98,102,221,48,253,226,152,37,
  179,16,145,34,136,54,208,148,206,143,150,219,189,241,210,19,92,131,
  56,70,64,30,66,182,163,195,72,126,110,107,58,40,84,250,133,186,61,202,
  94,155,159,10,21,121,43,78,212,229,172,115,243,167,87,7,112,192,247,
  140,128,99,13,103,74,222,237,49,197,254,24,227,165,153,119,38,184,180,
  124,17,68,146,217,35,32,137,46,55,63,209,91,149,188,207,205,144,135,151,
  178,220,252,190,97,242,86,211,171,20,42,93,158,132,60,57,83,71,109,65,
  162,31,45,67,216,183,123,164,118,196,23,73,236,127,12,111,246,108,161,59,
  82,41,157,85,170,251,96,134,177,187,204,62,90,203,89,95,176,156,169,160,
  81,11,245,22,235,122,117,44,215,79,174,213,233,230,231,173,232,116,214,
  244,234,168,80,88,175};
const int expq[255] = {1,2,4,8,16,32,64,128,29,58,116,232,205,135,19,38,76,
  152,45,90,180,117,234,201,143,3,6,12,24,48,96,192,157,39,78,156,
  37,74,148,53,106,212,181,119,238,193,159,35,70,140,5,10,20,40,80,
  160,93,186,105,210,185,111,222,161,95,190,97,194,153,47,94,188,101,
  202,137,15,30,60,120,240,253,231,211,187,107,214,177,127,254,
  225,223,163,91,182,113,226,217,175,67,134,17,34,68,136,13,26,52,104,
  208,189,103,206,129,31,62,124,248,237,199,147,59,118,236,197,151,51,
  102,204,133,23,46,92,184,109,218,169,79,158,33,66,132,21,42,84,168,
  77,154,41,82,164,85,170,73,146,57,114,228,213,183,115,230,209,191,99,
  198,145,63,126,252,229,215,179,123,246,241,255,227,219,171,75,150,49,
  98,196,149,55,110,220,165,87,174,65,130,25,50,100,200,141,7,14,28,56,
  112,224,221,167,83,166,81,162,89,178,121,242,249,239,195,155,43,86,172,
  69,138,9,18,36,72,144,61,122,244,245,247,243,251,235,203,139,11,22,44,
  88,176,125,250,233,207,131,27,54,108,216,173,71,142};
#endif


#if Q==2                        // please do not try this (i.e., Log2Q=1)
  #define GF_mul(a, b) ((a)&(b))
#else
int GF_mul(int a, int b)
{
  if (a == 0 || b == 0) return 0;
  if (a == 1) return b;
  if (b == 1) return a;
  return expq[(logq[a] + logq[b]) % (Q-1)];
}
#endif

#define GF_add(a, b) ((a)^(b))
#define GF_sub(a, b) ((a)^(b))

#ifdef LONGLONG
  typedef long long int NTT;    // 8-byte int (for both VC and gcc)
#else
  typedef int NTT;
#endif

// Followings are based on Prof. MacKay's FFT code, but do not work
// in my program, so I changed a little bit (see comments below)
void ntt(NTT p[Q])
{
  int b, factor = 1, rest;
  for (b = 0; b < Log2Q; b++) {
    for (rest = 0; rest < Q/2; rest++) {
      int restH = rest >> b;
      int restL = rest & (factor-1);
      int rest0 = (restH << (b+1)) + restL;
      NTT rest1 = rest0 + factor;
      NTT prest0 = p[rest0];
//      p[rest0] = p[rest1] - p[rest0];
//      p[rest1] += prest0;
      p[rest0] += p[rest1];
      p[rest1] = prest0 - p[rest1];
    }
    factor += factor;
  }
}

int n, m;
int rmax, cmax;
int *row_weight, *col_weight;
int **row_col, **row_N;
int **col_row, **col_N;
int **H;
int ***logra;
int ***logqa;
int ***logracol;
int ***logqacol;
NTT **fQa;
int **isnegative;
int *tmp_z;
int *tmp_x;

#define exp2(x) pow(2.0,x)
#define log2(x) (log(x)/log(2.0))

static unsigned int rndm = 2815UL;
void SRand(int n) {
  rndm = n;
//  srand(n);
}

#define RandMax 0x7fffffffUL
//#define RandMax RAND_MAX
unsigned int Rand(void)
{
//  return rand();
  return rndm = (77UL * rndm + 1243UL) & RandMax; // 31bit
}

#define INT 6/*8*/              // int part
#define DECI  14/*13*/              // fraction part
#define FMUL  (1<<DECI)       // multiplier
#define PREC  (1.0/FMUL)      // precision
#define LEVELS  (1<<(INT+DECI))
static int flog[LEVELS];
static unsigned int _fexp[LEVELS*2], *fexp = &_fexp[LEVELS];

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

#define exp(x) pow(2, x)
#define log(x) (log(x)/log(2.0))

void initlogexptab2(void)
{
  int i = 1;
  double right = log(fix2float(i) - 0.5*PREC);
  flog[0] = -FMUL*14;//-115000;
  for ( ; i < LEVELS; i++) {
    double d = fix2float(i);
    double left = log(d+0.5*PREC);
    flog[i] = float2fix((4*log(d)+right+left) / 6.0); // Simpson's law
    right = left;
  }
  i = -LEVELS;
  right = exp(fix2float(i) - 0.5*PREC);
  for ( ; i < LEVELS; i++) {
    double d = fix2float(i);
    double expd = exp(d);
    double left = exp(d+0.5*PREC);
    if (expd > (1UL<<(31-DECI)))
      fexp[i] = 1UL<<31;
    else
      fexp[i] = float2fixu((4*expd+right+left) / 6.0); // Simpson's law
    right = left;
  }
}

#if 0
unsigned int Fexp(int x)
{
  if (x < -LEVELS) return 0;
  else if (x >= LEVELS) {
    int i = 0;
    do {
      x -= FMUL;
      i++;
    } while (x >= LEVELS);
    return fexp[x] << i;
  }
  return fexp[x];
}
#else
  #define Fexp(x) fexp[x]
#endif

#if 1
int Flog(int x)
{
  int i;
  if (x <= 0) return -FMUL*14; //-115000
  for (i = 0; x >= LEVELS; i++)
    x >>= 1;
  return flog[x] + (i<<DECI);
}
#else
  #define Flog(x) flog[x]
#endif


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
  int i, j;
  for (j = 0; j < m; j++) {
    int sum = 0;
    for (i = 0; i < row_weight[j]; i++) {
      int xHmn = GF_mul(x[row_col[j][i]], H[j][i]);
      sum = GF_add(sum, xHmn);
    }
    s[j] = sum;
  }
}

NTT **malloc2DNTT(int a, int b) // allocates array[a][b]
{
  int i;
  NTT **pp = malloc(sizeof(NTT *) * a);
  NTT *p = malloc(sizeof(NTT) * a * b);
  if (pp == NULL || p == NULL) exit(-1);
  for (i = 0; i < a; i++) {
    pp[i] = p + b*i;
  }
  return pp;
}

int ***malloc2Dintp(int a, int b) // allocates array[a][b] (int pointer)
{
  int j;
  int ***ppp= malloc(sizeof(int **) * a);
  int **pp  = malloc(sizeof(int *) * a * b);
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

int ***malloc3Dint(int a, int b, int c) // allocates array[a][b][c]
{
  int i, j;
  int ***ppp= malloc(sizeof(int **) * a);
  int **pp  = malloc(sizeof(int *) * a * b);
  int *p    = malloc(sizeof(int) * a * b * c);
  if (ppp == NULL || pp == NULL || p == NULL) exit(-1);
  for (j = 0; j < a; j++) {
    for (i = 0; i < b; i++) {
      pp[i+b*j] = p + c*i + c*b*j;
    }
    ppp[j] = pp + b*j;
  }
  return ppp;
}

// Sum Product Decoder
// logfna: prior
// z: syndrome
// loop_max: max iteration
// x[]: original signal (just for reference)
int dec(int **logfna, int z[], int loop_max, int x[])
{
  int a, i, j, k, loop;
  int iir, prev = 999999, nodecr = 0;

  for (i = 0; i < n; i++) {
    for (k = 0; k < col_weight[i]; k++) {
      memcpy(logqacol[i][k], logfna[i], sizeof(int)*Q);
    }
  }

  for (loop = 0; loop < loop_max; loop++) {
    for (j = 0; j < m; j++) {
      int row_weightj = row_weight[j];
      int logprodqa[Q], sgnsum[Q];
      memset(sgnsum, 0, sizeof(sgnsum));
      memset(logprodqa, 0, sizeof(logprodqa));
      memset(*isnegative, 0, rmax * Q * sizeof(int));
      for (k = 0; k < row_weightj; k++) {
        for (a = 0; a < Q; a++) {
          fQa[k][GF_mul(a, H[j][k])] = Fexp(logqa[j][k][a]);
        }
        ntt(fQa[k]);
        //Flog(|fQa[k][a]|)をどこかに保存し後に使うのもアリ
        for (a = 0; a < Q; a++) {
          if (fQa[k][a] < 0) {
            fQa[k][a] = -fQa[k][a];
            isnegative[k][a] = 1;
            sgnsum[a] ^= 1;     // 負のときだけカウントアップ
          }
          fQa[k][a] = (fQa[k][a] + (1<<(Log2Q/2-1))) >> (Log2Q/2);
          logprodqa[a] += Flog(fQa[k][a]);
          assert(logprodqa[a] > -(1<<29));
        }
      }
      for (k = 0; k < row_weightj; k++) {
        NTT fRa[Q];
        for (a = 0; a < Q; a++) {
          fRa[a] = Fexp(logprodqa[a] - Flog(fQa[k][a]) + Log2Q*FMUL); //prod. except k
          if (isnegative[k][a] != sgnsum[a])
            fRa[a] = -fRa[a];
        }
        ntt(fRa);          // intt
        for (a = 0; a < Q; a++) {
          logra[j][k][a] = Flog(fRa[GF_sub(z[j], GF_mul(a, H[j][k]))]);
        }
      }
    }
    // update qa & tentative decode
    for (i = 0; i < n; i++) {
      int col = col_weight[i];
      int argmaxa = -1;
      int maxp = -(1<<30);
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
        int *qajk = logqa[j][k];
        int sum = 0;
        for (a = 0; a < Q; a++) {
          sum += Fexp(*qajk++);
        }
        if (sum <= 0) {//all zero happenes
          //printf("%d,%d:sum%d ", j, k, sum);
          sum = Q;//taka051207
        }
        sum = Flog(sum);
        qajk = logqa[j][k];
        for (a = 0; a < Q; a++) {
          *qajk++ -= sum;
        }
      }
    }
    //printf("%2d: HamDist(x,x^)=%d ", loop+1, HamDist(x, tmp_x, n));

    enc(tmp_x, tmp_z);
    {
      int dist = HamDist(z, tmp_z, m);
      //printf("HamDist(s,synd(x^))=%d\n", dist);
      if (dist == 0)           // nothing more can be done
        return 0;
      // nonconvergence detection
      if (loop == 0) iir = dist;
      else iir = (int)(iir * 0.85 + dist * 0.15 + 0.5);

      if (prev <= dist) nodecr++;
      else nodecr = 0;
      if (dist > iir * 1.1 || nodecr > 10) break; // no conversion
      prev = dist;
    }
  }

  return -1;
}

//http://www.inference.phy.cam.ac.uk/mackay/codes/alist.html
void initdec(char *s)
{
  int i, j, q, *count;
  int isbinary = 0;
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
    isbinary = 1;
    q = 2;
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
  fQa          = malloc2DNTT(rmax, Q);
  isnegative   = malloc2Dint(rmax, Q);

  {//skip n lines
    for (i = 0; i < n; i++) {
      for (j = 0; j < cmax; j++)
        fscanf(fp, isbinary ? "%*d" : "%*d%*d");
    }
  }

  for (j = 0; j < m; j++) {
    for (i = 0; i < row_weight[j]; i++) {
      int v;
      if (isbinary) {
        fscanf(fp, "%d", &v);
        H[j][i] = 1;
      } else
        fscanf(fp, "%d%d", &v, &H[j][i]);
      v--;
      row_col[j][i] = v;
      row_N[j][i] = count[v];
      col_row[v][count[v]] = j;
      col_N[v][count[v]] = i;
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

  logra  = malloc3Dint(m, rmax, Q);
  logqa  = malloc3Dint(m, rmax, Q);
  logracol=malloc2Dintp(n, cmax);
  logqacol=malloc2Dintp(n, cmax);
  for (i = 0; i < n; i++) {
    int k;
    int *cri = col_row[i], *cNi = col_N[i];
    int col = col_weight[i];
    for (k = 0; k < col; k++) {
      logracol[i][k] = logra[cri[k]][cNi[k]];
      logqacol[i][k] = logqa[cri[k]][cNi[k]];
    }
  }
}

int main(int argc, char **argv)
{
  int i, j, *s, *s0, *x, *y;
  int **logfna;
  double error;
  void (*channel)(int x[], int y[], double p, int **logfna) = lap;
  int iteration = 1000;

  initlogexptab2();

  if (argc == 6 && strcmp(argv[1], "-iter") == 0) {
    iteration = atoi(argv[2]);
    argc -= 2;
    argv += 2;
  }
  if (argc != 4) {
  usage:
    fprintf(stderr, "GFq_LDPC_NTT - GF(%d) simulator\n", Q);
    fprintf(stderr, "usage: GFq_LDPC_NTT [-iter num] <lap|bsc> NoiseLevel ParityMatrix\n");
    fprintf(stderr, "e.g.: GFq_LDPC_NTT lap 0.5 q8.sp.6000.4000.3000.1\n");
    fprintf(stderr, "      GFq_LDPC_NTT bsc 0.08 q4.sp.9000.6000.4500.1\n");
    return -1;
  }
  if (strcmp(argv[1], "bsc") == 0) channel = bsc;
  else if (strcmp(argv[1], "lap") != 0) {
    fprintf(stderr, "please specify 'lap' or 'bsc'\n");
    goto usage;
  }
  error = atof(argv[2]);
  initdec(argv[3]);
  s = malloc(sizeof(int) * m);  // syndrome
  s0= malloc(sizeof(int) * m);  // syndrome
  x = malloc(sizeof(int) * n);  // source
  y = malloc(sizeof(int) * n);  // side information
  logfna=malloc2Dint(n, Q);

  // Start Timer
  clock_t start = clock(), diff;

  // iterate experiments
  for (j = 1; j <= 3; j++) {
    //printf("\nGF(%d) experiment %d:\n", Q, j);
    SRand(j);
    for (i = 0; i < n; i++)
      x[i] = Rand() % Q;
    enc(x, s);
    channel(x, y, error, logfna);
    enc(y, s0);
    //printf("start: HamDist(x,y)=%d HamDist(synd(x),synd(y))=%d\n",
    //       HamDist(x,y,n), HamDist(s,s0,m));
    puts(
      dec(logfna, s, iteration, x) ?
      "failed." : "converged."
      );
  }

  // End Timer
  diff = clock() - start;
  int msec = diff * 1000 / CLOCKS_PER_SEC;
  printf("Time taken %d seconds %d milliseconds", msec/1000, msec%1000);

  return 0;
}
