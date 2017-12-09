/* GF(q>2) LDPC encoding/decoding simulator */
/* based on M. C. Davey et al. "Low-Density...over GF(q)" June 1998 */
/* 2005.11.23 by Seishi Takamura, Stanford University */
#include <mex.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#define Q 13
#define page_size 18336
#define QLC   4
#define CW_per_page   1
// If Q is defined, GF(Q) arithmetic is used.
// Else, let Log2Q defined and QF(2^Log2Q) arithmetic is used.

//#define Q 31  // must be prime, 31 works better than Q=32

#ifdef Q

  #define GF_mul(a, b)	(((a)*(b)) % Q)
  #define GF_add(a, b)	(((a)+(b)) % Q)
  #define GF_sub(a, b)	(((a)-(b)+Q*2) % Q)

#else // Q

#ifndef Log2Q
  #define Log2Q	5 // GF(2^5)
#endif
#define Q	(1<<Log2Q)                   // GF(Q)

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
#else
#error "illegal Q value"
#endif

int GF_mul(int a, int b)
{
  if (a == 0 || b == 0)  return 0;
  if (a == 1) return b;
  if (b == 1) return a;
  return expq[(logq[a] + logq[b]) % (Q-1)];
}

#define GF_add(a, b) ((a)^(b))
#define GF_sub(a, b) ((a)^(b))

#endif //Q

int grey_code_inv[16] = {6,5,7,14,9,12,8,13,3,4,2,15,10,11,1,0};

int n, m;
int rmax, cmax;
int *row_weight, *col_weight;
int **row_list_col, **row_list_N;
int **col_list_row, **col_list_N;
int **H;
double ***logra;
double ***logqa;
double ***logsigma;
double ***logrho;
int *tmp_z;
int *tmp_x;
double **p_sent_given_rec_T;


// Lower page data stored first
char data_read[QLC*page_size];
char data_written[QLC*page_size];


static unsigned int rndm = 2815UL;
void SRand(int n) {
//  rndm = n;
  srand(n);
}

//#define RandMax 0x7fffffffUL
#define RandMax RAND_MAX
unsigned int Rand(void)
{
    return rand();
//  return rndm = (77UL * rndm + 1243UL) & RandMax; // 31bit
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

void assign_llr(int y[], double* p_rec_given_sent, double **logfna, int num_reads){
  int i, j;

  for(i = 0; i < n; i++){
    for(j = 0; j < Q; j++){
      if(p_sent_given_rec_T[y[i]][j] == 0){
        logfna[j][i] = -100;
      } else{
        logfna[j][i] = log(p_sent_given_rec_T[y[i]][j]);
      }
    }
  }
}

void channel(int x[], int y[], double* p_rec_given_sent, double **logfna, int num_reads){
  int i, rand_select, rec_ind, j;
  double temp;
  
  for(i = 0; i < n; i++){
    temp = 0;
    rec_ind = 0;
    rand_select = rand()%101;

    if(!rand_select){
      while(!p_rec_given_sent[(rec_ind++)*Q+x[i]]);
    }

    while(temp < rand_select && rec_ind != Q*num_reads){
      temp += 100*p_rec_given_sent[(rec_ind++)*Q+x[i]];
    }

    y[i] = --rec_ind;

    for(j = 0; j < Q; j++){
      if(p_sent_given_rec_T[y[i]][j] == 0){
        logfna[j][i] = -100;
      } else{
        logfna[j][i] = log(p_sent_given_rec_T[y[i]][j]);
      }
    }
  }
}

// y[n] := x[n] + Laplacian_noise
// logfna[a][n]: log(Prob(x[n]=a | y[n]))
// stddev: noise level
void lap(int x[], int y[], double stddev, double **logfna)
{
  int i;
  int count[Q];
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
      double fna;
      //logfna[a][i] = -(abs(y[i] - a))/stddev;
      if (y[i] == a) logfna[a][i] = -0.5/stddev;//-0.2396/stddev; //l((1-e(-0.5))*2)
      else logfna[a][i] = -(abs(y[i] - a)-0.0413)/stddev; //l(e(0.5)-e(-0.5))
      fna = exp(logfna[a][i]);
      sum += fna;
    }
    //normalize
    sum = log(sum);
    for (a = 0; a < Q; a++) {
      logfna[a][i] -= sum;
    }
  }
  sum = 0;
  for (i = 0; i < Q; i++) {
    if (count[i]) sum += count[i] * log(count[i]);
  }
  printf("m/n=%g, ", (double)m/n);
  printf("noise entropy = %g bits, rate = %g\n", (-sum/n + log(n))/ log(2.0),
         (-sum/n + log(n)) / log(Q));
  printf("PSNR = %g\n", 10*log((Q-1)*(Q-1)/mse(x,y,n)));


/*  {int j;
    for (j = 0; j < n; j++) if (y[j] < Q-10) break;
    if (j < n) {
      FILE *fp = fopen("tmphist.dat", "wt");
      for (i = 0; i < 10; i++)
        fprintf(fp, "%d %d %g %g\n", i, count[i], logfna[y[j]+i][j],
		-i/stddev);
      fclose(fp);
    }
    exit(0);
  }*/
}

void enc(int x[], int s[])
{
  int i, j;
  for (j = 0; j < m; j++) {
    int sum = 0;
    for (i = 0; i < row_weight[j]; i++) {
      int xHmn = GF_mul(x[row_list_col[j][i]], H[j][i]);
      sum = GF_add(sum, xHmn);
    }
    s[j] = sum;
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
double ***malloc3Ddouble(int a, int b, int c) // allocates array[a][b][c]
{
  int i, j;
  double ***ppp= malloc(sizeof(double **) * a);
  double **pp  = malloc(sizeof(double *) * a * b);
  double *p    = malloc(sizeof(double) * a * b * c);
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
int dec(double **logfna, int z[], int loop_max)
{
  double sum;
  int a, i, j, k, s, t, loop;

  for (a = 0; a < Q; a++) {
    for (i = 0; i < n; i++) {
      for (k = 0; k < col_weight[i]; k++)
        logqa[a][col_list_row[i][k]][col_list_N[i][k]] = logfna[a][i];
    }
  }

  for (loop = 1; loop <= loop_max; loop++) {
    for (j = 0; j < m; j++) {
      // sigma
      for (k = 0; k < row_weight[j]; k++) {
        for (a = 0; a < Q; a++) {
          if (k == 0) {
            int aHmn = GF_mul(a, H[j][k]);
            logsigma[aHmn][j][k] = logqa[a][j][k];
          } else {
            sum = 0;
            for (t = 0; t < Q; t++) {
              s = GF_sub(a, GF_mul(t, H[j][k]));
              sum += exp(logsigma[s][j][k-1] + logqa[t][j][k]);
            }
            logsigma[a][j][k] = log(sum);
	        }
        }
      }
      // rho
      for (k = row_weight[j] - 1; k >= 0; k--) {
        if (k == row_weight[j]-1) {
          for (a = 0; a < Q; a++) {
            int aHmn = GF_mul(a, H[j][k]);
            logrho[aHmn][j][k] = logqa[a][j][k];
          }
        } else {
          for (a = 0; a < Q; a++) {
            sum = 0;
            for (t = 0; t < Q; t++) {
              s = GF_sub(a, GF_mul(t, H[j][k]));
              sum += exp(logrho[s][j][k+1] + logqa[t][j][k]);
            }
            logrho[a][j][k] = log(sum);
          }
        }
      }
    }
    // ra
    for (j = 0; j < m; j++) {
      for (k = 0; k < row_weight[j]; k++) {
        for (a = 0; a < Q; a++) {
          int z_aHmn = GF_sub(z[j], GF_mul(a, H[j][k]));
          if (k == 0) {
            logra[a][j][k] = logrho[z_aHmn][j][k+1];
          } else if (k == row_weight[j]-1) {
            logra[a][j][k] = logsigma[z_aHmn][j][k-1];
          } else {
            sum = 0;
            for (t = 0; t < Q; t++) {
              s = GF_sub(z_aHmn, t); // s = z[j] - a*Hmn - t
              sum += exp(logsigma[s][j][k-1] + logrho[t][j][k+1]);
            }
            logra[a][j][k] = log(sum);
          }
        }
      }
    }
    // update qa & tentative decode
    for (i = 0; i < n; i++) {
      int col = col_weight[i];
      int argmaxa = -1;
      double maxp = -9e99;
      for (a = 0; a < Q; a++) {
        double logprod = logfna[a][i];
        for (k = 0; k < col; k++) {
          logprod += logra[a][col_list_row[i][k]][col_list_N[i][k]];
        }
        for (k = 0; k < col; k++) {
          logqa[a][col_list_row[i][k]][col_list_N[i][k]] = logprod - logra[a][col_list_row[i][k]][col_list_N[i][k]]; // product except k
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
        sum = 0;
        for (a = 0; a < Q; a++) {
          sum += exp(logqa[a][j][k]);
        }
        sum = log(sum);
        for (a = 0; a < Q; a++) {
          logqa[a][j][k] -= sum;
        }
      }
    }
    
    enc(tmp_x, tmp_z);
    sum = HamDist(z, tmp_z, m);
    //printf("HamDist(s)=%g\n", sum);
    if (sum == 0)           // nothing more can be done
      return 0;
  }

  return 1;
}

void InitSumProductDecoderSyndrome(char *s)
{
  int i, j, q, *count;
  FILE *fp = fopen(s, "rt");
  if (fp == NULL) {
    fprintf(stderr, "cannot open %s\n", s);
    exit(-2);
  }
  fscanf(fp, "%d%d%d", &n, &m, &q);
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

  row_list_col = malloc2Dint(m, rmax);
  row_list_N   = malloc2Dint(m, rmax);
  col_list_row = malloc2Dint(n, cmax);
  col_list_N   = malloc2Dint(n, cmax);
  H            = malloc2Dint(m, rmax); // Hmn

  //skip n lines
  for (i = 0; i < n; i++) {
    for (j = 0; j < cmax; j++)
      fscanf(fp, "%*d%*d");
  }
  
  for (j = 0; j < m; j++) {
    for (i = 0; i < row_weight[j]; i++) {
      int v;
      fscanf(fp, "%d%d", &v, &H[j][i]);
      v--;
      row_list_col[j][i] = v;
      row_list_N[j][i] = count[v];
      col_list_row[v][count[v]] = j;
      col_list_N[v][count[v]] = i;
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
  // Assign random Hmn \in GF(Q)
  /*for (j = 0; j < m; j++) {
    for (i = 0; i < row_weight[j]; i++) {
      H[j][i] = 1 + (rand() % (Q-1)); // must prevent zero!!!
    }
  }*/
  
  free(count);
  fclose(fp);
  logra = malloc3Ddouble(Q, m, rmax);
  logqa = malloc3Ddouble(Q, m, rmax);
  logsigma=malloc3Ddouble(Q, m, rmax);
  logrho = malloc3Ddouble(Q, m, rmax);
  tmp_z = malloc(sizeof(int) * m);
  tmp_x = malloc(sizeof(int) * n);
}

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

// decode_mode: 1:  read data from file
//              0:  generate random data
void GFq_LDPC(int iteration, int num_trials, int num_reads, int decode_mode, double* p_rec_given_sent, double *errors)
{
  int i, j, k, l, block_sz_bytes, *s, *x, *y, dec_result;
  double **logfna;
  long data_ind;
  char lp_byte, mp_byte, up_byte, tp_byte;

  //InitSumProductDecoderSyndrome("bazar_GF_13.txt");
  n = 18336*8;
  block_sz_bytes = n/8;
  p_sent_given_rec_T = malloc2Ddouble(Q*num_reads, Q);

  s = malloc(sizeof(int) * m);  // syndrome
  x = malloc(sizeof(int) * n);  // source
  y = malloc(sizeof(int) * n);  // side information
  logfna=malloc2Ddouble(Q, n);

  /*printf("p_rec_given_sent = \n");
  for(i = 0; i < Q*num_reads; i++){
    for(j = 0; j < Q; j++){
      printf("%.2f ", p_rec_given_sent[i*Q+j]);
    }
    printf("\n");
  }*/

  make_p_sent_given_rec_T(p_rec_given_sent, num_reads);

  /*printf("p_sent_given_rec_T = \n");
  for(i = 0; i < Q*num_reads; i++){
    for(j = 0; j < Q; j++){
      printf("%.2f ", p_sent_given_rec_T[i][j]);
    }
    printf("\n");
  }*/

  errors[0] = 0;
  errors[1] = 0;

  // Start Timer
  clock_t start = clock(), diff;

  // iterate experiments
  if(decode_mode){
    //FILE *read_data = fopen("read_data.bin", "rb");
    FILE *written_data = fopen("written_data.bin", "rb");

    FILE *sym = fopen("sym.txt", "w");
    fprintf(sym, "%d ", 5);
    printf("hi\n");
    
    for(i = 0; i < num_trials/CW_per_page; i++){
      //fread(data_read, 1, QLC*page_size, read_data);
      fread(data_written, 1, QLC*page_size, written_data);

      for(j = 0; j < CW_per_page; j++){
        for(k = 0; k < block_sz_bytes; k++){
          lp_byte = data_written[j*block_sz_bytes+k];
          mp_byte = data_written[page_size+j*block_sz_bytes+k];
          up_byte = data_written[2*page_size+j*block_sz_bytes+k];
          tp_byte = data_written[3*page_size+j*block_sz_bytes+k];
          get_symbols_in_byte(x, k, lp_byte, mp_byte, up_byte, tp_byte);
          /*lp_byte = data_read[j*block_sz_bytes+k];
          mp_byte = data_read[page_size+j*block_sz_bytes+k];
          up_byte = data_read[2*page_size+j*block_sz_bytes+k];
          tp_byte = data_read[3*page_size+j*block_sz_bytes+k];
          get_symbols_in_byte(y, k, lp_byte, mp_byte, up_byte, tp_byte);*/
        }
        /*channel(x, y, p_rec_given_sent, logfna, num_reads);
        for(k = 0; k < n; k++){
          fprintf(sym, "%d ", y[k]);
        }*/

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
    
  } else{
    for (j = 1; j <= num_trials; j++) {
      srand(j);
      // Generate random data
      //printf("x = ");
      for (i = 0; i < n; i++){
        x[i] = rand() % Q;
      }

      enc(x, s);
      channel(x, y, p_rec_given_sent, logfna, num_reads);


      //lap(x, y, 0.01, logfna);
      dec_result = dec(logfna, s, iteration);
    
      if(dec_result){
        errors[0]++;
      } else {
        if(HamDist(tmp_x, x, n) != 0) errors[1]++;
      }
    }
  }
  

  // End Timer
  diff = clock() - start;
  int msec = diff * 1000 / CLOCKS_PER_SEC;
  printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);
}

// The gateway function
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *conf_mat;              // input scalar
    int num_trials;              // input scalar
    int max_iter;              // input scalar
    int num_reads;               // input scalar
    int decode_mode;          // input scalar
    
    // check for proper number of arguments
    if(nrhs!=5) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Five inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
    }
    
    // get the values of the inputs
    max_iter = mxGetScalar(prhs[0]);
    num_trials = mxGetScalar(prhs[1]);
    num_reads = mxGetScalar(prhs[2]);
    decode_mode = mxGetScalar(prhs[3]);
    conf_mat = mxGetPr(prhs[4]);
    
    // create the output matrix
    plhs[0] = mxCreateNumericMatrix(1,2,mxDOUBLE_CLASS,mxREAL);
    
    // get a pointer to the real data in the output matrix
    double *outMatrix = mxGetPr(plhs[0]);
    
    // call the computational routine
    GFq_LDPC(max_iter, num_trials, num_reads, decode_mode, conf_mat, outMatrix);
}
