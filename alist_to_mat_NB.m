% alist_to_mat
% Reads a NB LDPC matrix in alist format
% return a sparse representation of it

% Author: Mohammed Al Ai Baky
% Created: 12/9/2017

function H = alist_to_mat_NB(file)
    % Open the alist matrix file
    mat_file = fopen(file);
    
    % Create sparse matrix
    c_r_q = fscanf(mat_file, '%d', [1 3]);
    c = c_r_q(1);
    r = c_r_q(2);
    H = sparse(r,c);
    
    % Fetch the the max weights and weight vectors
    col_weight_max = fscanf(mat_file, '%d', 1);
    row_weight_max = fscanf(mat_file, '%d', 1);
    col_weight_arr = fscanf(mat_file, '%d', [1 c]);
    row_weight_arr = fscanf(mat_file, '%d', [1 r]);
    
    for i = 1:c
        col = fscanf(mat_file, '%d', [1 2*col_weight_arr(i)]);
        skip_size = 2*(col_weight_max - col_weight_arr(i));
        skip = fscanf(mat_file, '%d', [1 skip_size]);
        
        for j = 0:col_weight_arr(i)-1
            H(col(2*j+1),i) = col(2*j+2);
        end
    end
    
    fclose(mat_file);
    return;
end
