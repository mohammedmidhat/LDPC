function H = circulant_PEG(nvar, nchk, var_degree_sequence, chk_degree_sequence, p)
    H = sparse(nchk*p,nvar*p);
    
    cur_var_degree_sequence = zeros(1,p*nvar);
    cur_chk_degree_sequence = zeros(1,p*nvar);
    
    
end
