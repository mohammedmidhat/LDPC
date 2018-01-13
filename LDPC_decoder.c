#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <errno.h>
#include <time.h>

char hi_c[4];
void main(int argc, char **argv){
	int i;
	srand(time(NULL));
	for(i=0; i<10; i++){
		printf("%d ", rand()%100);
	}
}
