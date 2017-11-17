// GFq_LDPC_NTT.c, a GF(q>2) LDPC encoding/decoding simulator,
// powered by NTT (NumberTheoretic Transform, i.e., FFT over GF)
// Only allows GF(2^p), 2<=p<=8
// based on paper of M. C. Davey et al. "Low-Density...over GF(q)" June 1998
// (c) 2005-2006 by Seishi Takamura @ Stanford University / NTT (Nippon Telegraph and Telephone)
// Absolutely no warranty.
//#include <mex.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <errno.h>
#include <time.h>


#define Q   13
#define Log2Q 4
#define GF_add(a, b) (a+b)%Q
#define GF_sub(a, b) (a-b)%Q


// Inverse in GF(13), inv(a) = inv[a], inv[0] = 0 for consistency
int inv[13] = {0,1,7,9,10,8,11,2,5,3,4,6,12};

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
float **fQa, **fRa, **p_sent_given_rec_T;
int *tmp_z;
int *tmp_x;
int *conv_sequence;
float p_symbol_alloc;
float *p_symbol = &p_symbol_alloc;


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
      p_sent_given_rec_T[i][j] = P_x*(float) p_rec_given_sent[i][j]/P_y;
    }
  }
}

void find_all_config(int seq_sz_total, int seq_sz, int rest, int not_row){
  //printf("%d, %d, %d, %d\n", seq_sz_total, seq_sz, rest, not_row);
  int i, j;
  int acc;
  
  if(seq_sz == 1){
    int seq_iter = 0;
    float p_seq = 1.0;
    
    conv_sequence[seq_sz_total - seq_sz] = ((Q - rest)%Q+Q)%Q;
    //printf("seq = ");
    for(i = 0; i < seq_sz_total+1; i++){
      if(i != not_row){
        //printf("%d ", conv_sequence[seq_iter]);
        p_seq *= fQa[i][conv_sequence[seq_iter]];
        seq_iter++;
      }
    }
    *p_symbol += p_seq;
    //printf("\n");
    return;
  }
  for(i = 0; i < Q; i++){
    conv_sequence[seq_sz_total - seq_sz] = i;
    find_all_config(seq_sz_total, seq_sz-1, rest+i, not_row);
  }
}

void CNP(int row_weight, int syndrome){
  int i, a;
  
  for(i = 0; i < row_weight; i++){
    //printf("new config\n");
    for (a = 0; a < Q; a++){
      //printf("config[%d]\n", a);
      *p_symbol = 0;
      find_all_config(row_weight-1, row_weight-1, (((a-syndrome)%Q)+Q)%Q, i);
      fRa[i][a] = *p_symbol;
    }
  }
}

void channel(int x[], int y[], double** p_rec_given_sent, float **logfna){
  int i, rand_select, rec_ind, j;
  float temp;

  for(i = 0; i < n; i++){
    temp = 0;
    rec_ind = 0;
    rand_select = rand()%101;
    //printf("rand_select = %d\n", rand_select);
    while(temp <= rand_select && rec_ind != 13){
      //printf("temp = %.2f\n", temp);
      temp += 100*(float) p_rec_given_sent[rec_ind++][x[i]];
    }
    y[i] = --rec_ind;
    logfna[i] = p_sent_given_rec_T[y[i]];
    
    /*printf("y[%d] = %d\n", i, rec_ind);
    for(j = 0; j < Q; j++){
      printf("%.2f ", logfna[i][j]);
    }
    printf("\n");*/
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

// a: number of rows
// b: number of columns
void set_2D_float_arr_one(float **arr, int a, int b){
  int i, j;
  for(i = 0; i < a; i++){
    for(j = 0; j < b; j++){
      arr[i][j] = 1.0;
    }
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

double** malloc_2D_double(int m, int n){
  int i;
  double **pp = malloc(sizeof(double *) * m);
  double *p = malloc(sizeof(double) * m * n);
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
int dec(float **logfna, int z[], int loop_max)
{
  int a, i, j, k, l, loop;
  float normalizing_val;
  //printf("exp = %.2f\n",logfna[26]);
  //memcpy(logqacol[27][0], logfna[27], sizeof(float)*Q);
  //printf("p-vec = \n");
  for (i = 0; i < n; i++) {
    /*for(j = 0; j < Q; j++){
      printf("%.2f ", logfna[i][j]);
    }
    printf("\n");*/
    for (k = 0; k < col_weight[i]; k++) {
      memcpy(logqacol[i][k], logfna[i], sizeof(float)*Q);
    }
  }
  //printf("hi\n");

  for (loop = 0; loop < loop_max; loop++) {
    for (j = 0; j < m; j++) {
      int row_weightj = row_weight[j];
      float logprodqa[Q];
      //int sgnsum[Q];
      //memset(sgnsum, 0, sizeof(sgnsum));
      memset(logprodqa, 1.0, sizeof(logprodqa));
      set_2D_float_arr_one(fRa, rmax, Q);
      /*if(j < 1 && loop == 0){
        printf("check = %d\n", j);
        printf("row_weightj = %d\n", row_weightj);
      }*/
      for (k = 0; k < row_weightj; k++) {
        for (a = 0; a < Q; a++) {
          fQa[k][GF_mul(a, H[j][k])] = logqa[j][k][a];
        }
        
        /*if(j < 1 && loop == 0){
          printf("fQa[%d] = ", k);
          for(a = 0; a < Q; a++){
            printf("%.2f,", fQa[k][a]);
          }
          printf("\n");
        }*/
      }

      CNP(row_weightj, z[j]);

      for (k = 0; k < row_weightj; k++) {
        normalizing_val = 0;
        /*if(j < 1 && loop == 0){
          printf("logra[%d][%d] = ", j, k);
        }*/
        for (a = 0; a < Q; a++) {
          logra[j][k][a] = fRa[k][GF_mul(a, H[j][k])];
          normalizing_val += logra[j][k][a];
          /*if(j < 1 && loop == 0){
            printf("%.2f,", logra[j][k][a]);
          }*/
        }
        /*if(j < 1 && loop == 0){
          printf("\n");
        }*/
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
        float logprod = logfna[i][a];
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
        }
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
  fQa          = malloc_2D_float(rmax, Q);
  fRa          = malloc_2D_float(rmax, Q);
  conv_sequence = malloc((rmax-1)*sizeof(int));

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
/*
void test_NB_LDPC_p(int iteration, int num_trials, double** p_rec_given_sent, double *errors)
{
  printf('hi');
  int i, j, *s, *x, *y, dec_result;
  float **logfna;

  printf('hi');
  initdec("bazar_GF_13.txt");
  printf('hi');
  p_sent_given_rec_T = malloc_2D_float(Q, Q);
  printf('hi');
  s = malloc(sizeof(int) * m);  // syndrome
  x = malloc(sizeof(int) * n);  // source
  y = malloc(sizeof(int) * n);  // side information
  logfna = malloc_2D_float(n, Q);
  printf('hi');
  make_p_sent_given_rec_T(p_rec_given_sent);
  printf('hi');
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
*/
int main(int argc, char **argv)
{
  int i, j, *s, *x, *y, dec_result;
  float **logfna;


  int iteration = 30;
  int num_trials = 1;
  double errors[2];
  double** p_rec_given_sent;
  p_rec_given_sent = malloc_2D_double(Q, Q);
  memset(p_rec_given_sent[0], 0, 13*sizeof(double));
  memset(p_rec_given_sent[1], 0, 13*sizeof(double));
  memset(p_rec_given_sent[2], 0, 13*sizeof(double));
  memset(p_rec_given_sent[3], 0, 13*sizeof(double));
  memset(p_rec_given_sent[4], 0, 13*sizeof(double));
  memset(p_rec_given_sent[5], 0, 13*sizeof(double));
  memset(p_rec_given_sent[6], 0, 13*sizeof(double));
  memset(p_rec_given_sent[7], 0, 13*sizeof(double));
  memset(p_rec_given_sent[8], 0, 13*sizeof(double));
  memset(p_rec_given_sent[9], 0, 13*sizeof(double));
  memset(p_rec_given_sent[10], 0, 13*sizeof(double));
  memset(p_rec_given_sent[11], 0, 13*sizeof(double));
  memset(p_rec_given_sent[12], 0, 13*sizeof(double));
  //printf("hi\n");
  for(i = 0; i < Q; i++){
    p_rec_given_sent[i][i] = 1.00;
  }
  /*for(i = 1; i < Q; i++){
    p_rec_given_sent[i][i-1] = 0.02;
  }
  p_rec_given_sent[11][12] = 0.02;*/



  initdec("bazar_GF_13.txt");
  p_sent_given_rec_T = malloc_2D_float(Q, Q);
  
  s = malloc(sizeof(int) * m);  // syndrome
  x = malloc(sizeof(int) * n);  // source
  y = malloc(sizeof(int) * n);  // side information
  logfna = malloc_2D_float(n, Q);

  make_p_sent_given_rec_T(p_rec_given_sent);
  /*printf("p_sent_given_rec_T = \n");
  for(i = 0; i < Q; i++){
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
  for (j = 1; j <= num_trials; j++) {
    srand(j);
    // Generate random data
    for (i = 0; i < n; i++){
      x[i] = rand() % Q;
    }
    
    enc(x, s);
    channel(x, y, p_rec_given_sent, logfna);
    /*printf("y = ");
    for(i = 0; i < n; i++){
      printf("%d ", y[i]);
    }
    printf("\n s = ");
    for(i = 0; i < m; i++){
      printf("%d ", s[i]);
    }
    printf("\n");*/
    dec_result = dec(logfna, s, iteration);
    
    if(dec_result){
      errors[0]++;
    } else {
      if(HamDist(tmp_x, x, n) != 0) errors[1]++;
    }
  }

  printf("errors = %.2f\n", errors[0]);
  printf("undet_errors = %.2f\n", errors[1]);

  // End Timer
  diff = clock() - start;
  int msec = diff * 1000 / CLOCKS_PER_SEC;
  printf("Time taken %d seconds %d milliseconds", msec/1000, msec%1000);

  return 0;
}
/*
// The gateway function
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double **conf_mat;              // input scalar
    int num_trials;              // input scalar
    int max_iter;              // input scalar
    
    // check for proper number of arguments
    if(nrhs!=3) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Three inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
    }
    
    // get the values of the inputs
    max_iter = mxGetScalar(prhs[0]);
    num_trials = mxGetScalar(prhs[1]);
    conf_mat = mxGetPr(prhs[2]);
    
    // create the output matrix
    plhs[0] = mxCreateNumericMatrix(1,2,mxDOUBLE_CLASS,mxREAL);
    
    // get a pointer to the real data in the output matrix
    double *outMatrix = mxGetPr(plhs[0]);
    printf('hi');
    // call the computational routine
    test_NB_LDPC_p(max_iter, num_trials, conf_mat, outMatrix);
}*/
