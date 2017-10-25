function H = parse_MacKay_NB(file, k, n)
    % Read matrix text file
    array = dlmread(file);
    
    % Initializing a sparse matrix according to the matrix text file read
    H = spalloc(n-k, n, nnz(array)/2);
    
    % Getting the row size
    [dummy,row_size] = size(array(1,:));
    
    % Fill in the sparse matrix according to the matrix text file read
    for vr = 1:n
        for cn = 1:2:row_size
            if(array(vr, cn) ~= 0)
                H(array(vr,cn),vr) = array(vr,cn+1);
            end
        end
    end
    
    return;
end
