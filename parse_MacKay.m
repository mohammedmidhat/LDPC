% parse_MacKay
% Parse the parity-check matrix (H) in the alist format into a sparse
% matrix

% Author: Mohammed Al Ai Baky
% Created: 8/10/2017

function H = parse_MacKay(file, k, n)
    array = dlmread(file);
    H = spalloc(n-k, n, prod(size(array)));
    for i = 1:n
        H(array(i, :),i) = 1;
    end
    return;
end
