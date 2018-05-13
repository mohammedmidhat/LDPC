% sort_alist
% Sort the indices in the alist file of the LDPC code

% Author: Mohammed Al Ai Baky
% Created: 2/22/2018


function sort_alist(file_name_r, file_name_w)
    file_ID_r = fopen(file_name_r, 'r');
    file_ID_w = fopen(file_name_w, 'w');
    
    c_r = fscanf(file_ID_r, '%d', [1 2]);
    fprintf(file_ID_w,'%d ', c_r);
    fprintf(file_ID_w,'\n');
    
    col_row_weights_max = fscanf(file_ID_r, '%d', [1 2]);
    fprintf(file_ID_w,'%d ', col_row_weights_max);
    fprintf(file_ID_w,'\n');
    
    col_weight_arr = fscanf(file_ID_r, '%d', [1 c_r(1)]);
    row_weight_arr = fscanf(file_ID_r, '%d', [1 c_r(2)]);
    fprintf(file_ID_w,'%d ', col_weight_arr);
    fprintf(file_ID_w,'\n');
    fprintf(file_ID_w,'%d ', row_weight_arr);
    fprintf(file_ID_w,'\n');
    
    for i = 1:c_r(1)
        indices = fscanf(file_ID_r, '%d', [1 col_weight_arr(i)]);
        indices = sort(indices);
        fprintf(file_ID_w,'%d ', indices);
        fprintf(file_ID_w,'%d ', zeros(1,col_row_weights_max(1)-col_weight_arr(i)));
        fprintf(file_ID_w,'\n');
    end
    for i = 1:c_r(2)
        indices = fscanf(file_ID_r, '%d', [1 row_weight_arr(i)]);
        indices = sort(indices);
        fprintf(file_ID_w,'%d ', indices);
        fprintf(file_ID_w,'%d ', zeros(1,col_row_weights_max(2)-row_weight_arr(i)));
        fprintf(file_ID_w,'\n');
    end
    
    fclose(file_ID_w);
    fclose(file_ID_r);
end
