#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <errno.h>
#include <time.h>


void main(int argc, char **argv){
	int x, y;
	x = 220;
	y = 260;

    FILE *fil = fopen("exp.bin", "wb");

    fwrite(&x, 1, 1, fil);
    fwrite(&y, 1, 1, fil);
}
