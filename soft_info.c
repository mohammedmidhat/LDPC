#define page_size	18336
#define QLC		4


FILE *file_pointer_arr;



void main(void){
	int i, j;
	char string_data[page_size*QLC*NUM_R];

	file_pointer_arr = malloc(NUM_R*sizeof(FILE *));
	for(i = 0; i < NUM_R; i++){
		file_pointer_arr[i] = fopen(argv[i+1], "rb");
	}

	for(i = 0; i < NUM_R; i++){
		for(j = 0; j < QLC; j++){
			fread(string_data+page_size*(i*QLC+j), page_size, sizeof(char), file_pointer_arr[i]);
		}
	}

	
}
