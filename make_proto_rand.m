function H = make_proto_rand(n,k,H_proto)
    [num_circ_row,num_circ_col] = size(H_proto);
    circ_size = n/num_circ_col;
    H = spalloc(n-k,n,circ_size*nnz(H_proto));
    
    for i = 1:num_circ_col
        for l = 1:num_circ_row
            permutation = randperm(circ_size);
            for j = 1:circ_size
                col_offset = (j-1)*num_circ_col;
                H((permutation(j)-1)*num_circ_row+l,i+col_offset) = 1;
            end
        end
    end                
end
