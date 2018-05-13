/*
Creates LDPC code with PEG method (Not quasi-cyclic)
nvar and nchk defines the number of variable and check nodes respectively
degree_sequence is an array of the variable node degrees

The command line argument passed in is text file to store the H-matrix in the alist format

Author: Mohammed Al Ai Baky
Created: 9/24/2017
*/


#include <queue>          // std::queue
#include <stdio.h>
#include <stdlib.h>


int nvar, nchk;
int *degree_sequence, *chk_degrees;
int *H;


int count_ones(int array[], int size){
	int i, num_ones = 0;
	for(i = 0; i < size; i++){
		num_ones += array[i];
	}
	return num_ones;
}

int array_equal(int array_1[], int array_2[], int size){
	int i;
	for(i = 0; i < size; i++){
		if(array_1[i] != array_2[i]){
			return 0;
		}
	}
	return 1;
}

int find_index_smallest(int array[], int low, int high){
	int mid, min_left, min_right;

	if(low == high){
		return 0;
	}
	if(high == low + 1){
		if (array[low] > array[high]){
			return high;
		} else{
			return low;
		}
	}
	
	mid = (low + high)/2;
	min_left = find_index_smallest(array, low, mid);
	min_right = find_index_smallest(array, mid+1, high);
	if(array[min_left] > array[min_right]){
		return min_right;
	} else{
		return min_left;
	}
}

int find_max_array(int array[], int low, int high){
	int mid, min_left, min_right;

	if(low == high){
		return array[low];
	}
	if(high == low + 1){
		if (array[low] > array[high]){
			return array[low];
		} else{
			return array[high];
		}
	}
	
	mid = (low + high)/2;
	min_left = find_max_array(array, low, mid);
	min_right = find_max_array(array, mid+1, high);
	if(min_left > min_right){
		return min_left;
	} else{
		return min_right;
	}
}

int find_smallest_chk(int array[], int size){
	int i, index_cur_chk_list_comp, size_cur_chk_list_comp = 0;

	int index[nchk];
	int degree[nchk];

	for(i = 0; i < size; i++){
		if(array[i] == 0){
			index[size_cur_chk_list_comp] = i;
			degree[size_cur_chk_list_comp] = chk_degrees[i];
			size_cur_chk_list_comp++;
		}
	}
	index_cur_chk_list_comp = find_index_smallest(degree, 0, size_cur_chk_list_comp-1);
	return index[index_cur_chk_list_comp];
}

int bfs(int var){
	int i, new_chk;

	int var_list[nvar] = {0};
	var_list[var] = 1;
	int cur_chk_list[nchk] = {0};
	int new_chk_list[nchk] = {0};
	std::queue<int> chk_Q;
	std::queue<int> var_Q;
	var_Q.push(var);

	while(1){
		while(!var_Q.empty()){
			for(i = 0; i < nchk; i++){
				if(H[i*nvar+var_Q.front()]){
					if(cur_chk_list[i]==0){
						new_chk_list[i] = 1;
						chk_Q.push(i);
					}
				}
			}
			var_Q.pop();
		}
		while(!chk_Q.empty()){
			for(i = 0; i < nvar; i++){
				if(H[chk_Q.front()*nvar+i]){
					if(var_list[i]==0){
						var_list[i] = 1;
						var_Q.push(i);
					}
				}
			}
			chk_Q.pop();
		}
		if(count_ones(new_chk_list, nchk) == nchk){
			new_chk = find_smallest_chk(cur_chk_list, nchk);
			return new_chk;
		} else if(array_equal(new_chk_list, cur_chk_list, nchk)){
			new_chk = find_smallest_chk(cur_chk_list, nchk);
			return new_chk;
		} else{
			for(i = 0; i < nchk; i++){
				cur_chk_list[i] = new_chk_list[i];
			}
		}
	}
}


int main(int argc, char **argv){
	int i, j, k, var, smallest_degree_chk, chk;
	int col_weight_max, row_weight_max, row_weight_counter;

	nvar = 204;
	nchk = 102;

	degree_sequence = (int*) malloc(nvar*sizeof(int));
	for(i = 0; i < nvar; i++){
		degree_sequence[i] = 3;
	}

	chk_degrees = (int*) malloc(nchk*sizeof(int));
	for(i = 0; i < nchk; i++){
		chk_degrees[i] = 0;
	}

	H = (int*) malloc(nvar*nchk*sizeof(int));
	for(i = 0; i < nchk; i++){
		for(j = 0; j < nvar; j++){
			H[i*nvar+j] = 0;
		}
	}


	// progressive_edge_growth
	for(var = 0; var < nvar; var++){
		printf("edge growth at var %d\n", var);
		for(k = 0; k < degree_sequence[var]; k++){
			if(k == 0){
				smallest_degree_chk = find_index_smallest(chk_degrees, 0, nchk-1);
				H[smallest_degree_chk*nvar+var] = 1;
				chk_degrees[smallest_degree_chk] += 1;
			} else{
				chk = bfs(var);
				H[chk*nvar+var] = 1;
				chk_degrees[chk] += 1;
			}
		}
	}


	// Write H into a .txt file in alist format
	FILE *data_sym_file = fopen(argv[1], "w");

	fprintf(data_sym_file, "%d %d\n", nvar, nchk);

	col_weight_max = degree_sequence[0];
	row_weight_max = find_max_array(chk_degrees , 0, nchk-1);
	fprintf(data_sym_file, "%d %d\n", col_weight_max, row_weight_max);
	
	for(i = 0; i < nvar; i++){
		fprintf(data_sym_file, "%d ", degree_sequence[i]);
	}
	
	fprintf(data_sym_file, "\n");
	for(i = 0; i < nchk; i++){
		fprintf(data_sym_file, "%d ", chk_degrees[i]);
	}
	fprintf(data_sym_file, "\n");
	
	for(i = 0; i < nvar; i++){
		for(j = 0; j < nchk; j++){
			if(H[j*nvar+i]){
				fprintf(data_sym_file, "%d ", j+1);
			}
		}
		fprintf(data_sym_file, "\n");
	}
	
	for(i = 0; i < nchk; i++){
		row_weight_counter = 0;
		for(j = 0; j < nvar; j++){
			if(H[i*nvar+j]){
				fprintf(data_sym_file, "%d ", j+1);
				row_weight_counter++;
			}
		}
		for(; row_weight_counter < row_weight_max; row_weight_counter++){
			fprintf(data_sym_file, "%d ", 0);
		}
		fprintf(data_sym_file, "\n");
	}
	fclose(data_sym_file);

	return 0;
}
