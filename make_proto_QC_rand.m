% make_proto_QC_rand
% Creates a binary LDPC parity-check matrix (H)
% The H-matrix is Quasi-Cyclic Protograph-based with random permutations

% Issues
% 1- The permutations must be modified to result in an actual Quasi-Cyclic
% matrix

% Author: Mohammed Al Ai Baky
% Created: 9/28/2017

function H = make_proto_QC_rand(n,k,H_proto)
    [num_circ_row,num_circ_col] = size(H_proto);
    circ_size = n/num_circ_col;
    H = spalloc(n-k,n,circ_size*nnz(H_proto));
    
    for j = 1:num_circ_col
        col_offset = (j-1)*circ_size;
        for i = 1:num_circ_row
            if(H_proto(i,j))
                row_offset = (i-1)*circ_size;
                permutation = randperm(circ_size);
                for k = 1:circ_size
                    H(permutation(1,k)+row_offset,k+col_offset) = 1;
                end
            end
        end
    end                
end
