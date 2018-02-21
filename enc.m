% enc
% Compute the ECC parity of of a codeword
% Information passed to it from initdec()
% Implements Takamura's encoding algorithm

% Author: Mohammed Al Ai Baky
% Created: 2/17/2018

function s = enc(x, row_weight, m, row_list_col)
    s = zeros(1,m);
    
    for j = 1:m
        for i = 1:row_weight(1,j)
            s(1,j) = xor(s(1,j),x(1,row_list_col(j,i)));
        end
    end

    return;
end
