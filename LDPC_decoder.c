#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <errno.h>
#include <time.h>

 
 void printPar(int l, int r, char* str, int count){
 	if(l < 0 || r < l) return;
 	if(l == 0 && r == 0){
 		printf("%s ", str);
 	} else{
 		if(l > 0){
 			*(str+count) = (char) '(';
 			printPar(l-1, r, str, count+1);
 		}
 		if(r>l){
 			*(str+count) = (char) ')';
 			printPar(l, r-1, str, count+1);
 		}
 	}
 }


void printP(int count){
	char str[6];
	printPar(count, count, str, 0);
}

int main(void){
	//printP(3);
	
	return 0;
}
