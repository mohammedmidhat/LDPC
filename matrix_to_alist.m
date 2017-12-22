% matrix_to_alist
% Writes a sparse matrix in the alist format into a text file
% The alist format is as in Prof. MacKay's ldpc code library

% Author: Mohammed Al Ai Baky
% Created: 9/24/2017


function matrix_to_alist(H, file_name)
    file_ID = fopen(file_name, 'w');
    
    [r,c] = size(H);
    fprintf(file_ID,'%d %d\n', c, r);
    
    col_weights_arr = zeros(1,c);
    row_weights_arr = zeros(1,r);
    for i = 1:c
        col_weights_arr(1,i) = nnz(H(:,i));
    end
    for i = 1:r
        row_weights_arr(1,i) = nnz(H(i,:));
    end
    col_weight_max = max(col_weights_arr);
    row_weight_max = max(row_weights_arr);
    fprintf(file_ID,'%d %d\n', col_weight_max, row_weight_max);
    fprintf(file_ID,'%d ', col_weights_arr);
    fprintf(file_ID,'\n');
    fprintf(file_ID,'%d ', row_weights_arr);
    fprintf(file_ID,'\n');
    
    for i = 1:c
        fprintf(file_ID,'%d ',find(H(:,i)));
        for j = 1:col_weight_max - nnz(find(H(:,i)))
            fprintf(file_ID,'%d ', 0);
        end
        fprintf(file_ID,'\n');
    end
    for i = 1:r
        fprintf(file_ID,'%d ',find(H(i,:)));
        for j = 1:row_weight_max - nnz(find(H(i,:)))
            fprintf(file_ID,'%d ', 0);
        end
        fprintf(file_ID,'\n');
    end
    
    fclose(file_ID);
    
end
