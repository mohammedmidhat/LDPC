% InitSumProductDecoderSyndrome
% Parses the H-matrix in the alist format
% and passes the information through the GFq_LDPC() mex call

% Author: Mohammed Al Ai Baky
% Created: 1/28/2018

function [n, m, cmax, rmax, col_weight, row_weight, row_list_col, col_list_row, col_list_N, H] = InitSumProductDecoderSyndrome(filename)
    code_file = fopen(filename);
    
    n_m_q = fscanf(code_file, '%d', [1 3]);
    n = n_m_q(1);
    m = n_m_q(2);
    q = n_m_q(3);
    cmax_rmax = fscanf(code_file, '%d', [1 2]);
    cmax = cmax_rmax(1);
    rmax = cmax_rmax(2);
    col_weight = fscanf(code_file, '%d', [1 n]);
    row_weight = fscanf(code_file, '%d', [1 m]);
    
    count = zeros(1,n);
    
    row_list_col = zeros(m, rmax);
    row_list_N = zeros(m, rmax);
    col_list_row = zeros(n, cmax);
    col_list_N = zeros(n, cmax);
    H = zeros(m, rmax);
    
    %skip n lines
    for i = 1:n
        skip = fscanf(code_file, '%d', [1 2*cmax]);
    end
    
    for j = 1:m
        for i = 1:row_weight(1,j)
            v = fscanf(code_file, '%d', 1);
            H(j, i) = fscanf(code_file, '%d', 1);
            
            row_list_col(j,i) = v-1;      %
            row_list_N(j,i) = count(1,v);
            col_list_row(v,count(1,v)+1) = j-1;
            col_list_N(v,count(1,v)+1) = i-1;
            count(1,v) = count(1,v) + 1;
        end
        
        skip_size = rmax - row_weight(1,j);
        skip = fscanf(code_file, '%d', [1 2*skip_size]);
    end
    
    fclose(code_file);
    return;
end
