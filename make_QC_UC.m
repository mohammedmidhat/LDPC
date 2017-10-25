H = parse_MacKay('204.33.484 (N=204,K=102,M=102,R=0.5).txt', 102, 204);

[num_rows,num_cols] = size(H);
[row_ind,col_ind] = find(H);
vr_deg = nnz(H(:,1));

for col = 1:num_rows/vr_deg
    start_ind = vr_deg*(col-1)+1;
    end_ind = vr_deg*col;
    candidate_rows = row_ind((start_ind:end_ind),1);
    H((col,),:) = H();