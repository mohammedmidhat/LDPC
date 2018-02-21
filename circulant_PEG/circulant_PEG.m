% circulant_PEG
% Generates LDPC code by circulant PEG (Quasi-cyclic)

% Author: Mohammed Al Ai Baky
% Created: 1/11/2018

function H = circulant_PEG(nvar, nchk, var_degree_sequence, chk_degree_sequence, p)
    H = sparse(nchk*p,nvar*p);
    
    cur_var_degree_sequence = zeros(1,p*nvar);
    cur_chk_degree_sequence = zeros(1,p*nvar);
    
    
end
