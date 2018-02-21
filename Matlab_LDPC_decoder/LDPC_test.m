H = alist_to_mat('204.33.484.txt');
[m n] = size(H);
bit_err_prob = [0.007:-0.0005:0.005];
err_count_array = zeros(1,size(bit_err_prob,2));
undet_err_array = zeros(1,size(bit_err_prob,2));
num_trials = 10;
tic;
times = zeros(size(bit_err_prob, 2), num_trials);
for j = 1:size(bit_err_prob,2)
    err_count = 0;
    undet_err = 0;
    for (k = 1:num_trials)
        codeword = rand(1, n)<0.5;
        noise = rand(1, size(codeword,2))<bit_err_prob(1,j);
        codeword_rec = mod(codeword+noise,2);
        
        llr = zeros(1,size(codeword_rec,2));
        llr(codeword_rec == 0) = log((1-bit_err_prob(1,j))/bit_err_prob(1,j));
        llr(codeword_rec == 1) = log(bit_err_prob(1,j)/(1-bit_err_prob(1,j)));
        syndrome = mod(H*codeword',2);
        
        [cw num_iter suc] = mpdec(H, llr, 30, syndrome);
        if(~isequal(codeword,cw))
            err_count = err_count + 1;
            if(suc)
                undet_err = undet_err + 1;
            end
        end
    end
    err_count_array(1,j) = err_count;
    undet_err_array(1,j) = undet_err;
end
t_post = toc;
loglog(bit_err_prob, err_count_array);
hold on;
loglog(bit_err_prob, undet_err_array);
title('(204,102) LDPC, 30 trials, min-sum, max iterations = 20');
xlabel('BER');
ylabel('Error Count');
set(gca, 'xdir', 'reverse');