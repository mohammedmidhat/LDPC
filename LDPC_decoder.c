#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <errno.h>
#include <time.h>


void main(void){
    /*FILE *fp = fopen("sn.bin", "rb");
    char buf[3];
    fgets(buf, 3, fp);

    for(int i = 0; i < 3; i++){
        printf("%d ", buf[i]);
    }*/
    int i;
    float arr[4];
    //memset(arr, 1.0, 4);
    for(i = 0; i < 4; i++){
        arr[i] = 1.0;
    }
    for(i = 0; i < 4; i++){
        printf("%.2f ", arr[i]);
    }
}
