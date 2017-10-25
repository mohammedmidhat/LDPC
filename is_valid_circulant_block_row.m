function [success, row_indices, col_indices] = is_valid_circulant_block_row(mat)
    
    success = 1;
    row_indices = [];
    col_indices = [];
    
    [r,c] = size(mat);
    last_row = mat(r,:);
    cn_deg = nnz(mat(r,:));
    
    for i = 1:r-1
        if nnz(mod(last_row+mat(i,:),2)) ~= 2*cn_deg;
            row_indices = [row_indices i];
            col_indices = [col_indices find(find(last_row)==find(mat(i,:)))];
        else
            return;
        end
    end
    
    success = 0;    
    return;
    
end % fxn