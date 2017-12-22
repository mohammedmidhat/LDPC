num_trials = 256;
max_num_iter = 30;
p = [0.24:-0.01:0.2];
err_count = zeros(1,size(p,2));
undet_err = zeros(1,size(p,2));

tic;
parfor (j = 1:size(p,2),4)
    result = GFq_LDPC_lap(max_num_iter,num_trials,p(1,j));
    err_count(1,j) = result(1)/num_trials;
    undet_err(1,j) = result(2)/num_trials;
end
t_post = toc;

loglog(p, err_count);
hold on;
loglog(p, undet_err,'--o');
