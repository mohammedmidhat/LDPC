%H = parse_MacKay('204.33.484 (N=204,K=102,M=102,R=0.5).txt', 102, 204);
G_t = null2(H);
%bit_err_prob = [0.05];
bit_err_prob = [0.1:-0.01:0.01];
err_count_array = zeros(1,size(bit_err_prob,2));
undet_err_array = zeros(1,size(bit_err_prob,2));
num_trials = 1000;
tic;
times = zeros(size(bit_err_prob, 2), num_trials);
for j = 1:size(bit_err_prob,2)
    err_count = 0;
    undet_err = 0;
    parfor (k = 1:num_trials,4)
        codeword = mod(sum(G_t(:,rand(1, size(G_t,2))<0.5),2),2)';
        noise = rand(1, size(codeword,2))<bit_err_prob(1,j);
        %noise = zeros(1,size(codeword,2));
        %noise(1,1:2) = 1;
        codeword_rec = mod(codeword+noise,2);
        %codeword_rec = noise;
        llr = zeros(1,size(codeword_rec,2));
        llr(codeword_rec == 0) = log((1-bit_err_prob(1,j))/bit_err_prob(1,j));
        llr(codeword_rec == 1) = log(bit_err_prob(1,j)/(1-bit_err_prob(1,j)));
        
        %t1 = toc;
        [cw num_iter, suc] = mpdec(H, llr, 30);
        %times(j, k) = toc - t1;
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
figure;
hold on;
plot(bit_err_prob, err_count_array);
plot(bit_err_prob, undet_err_array);
title('(204,102) LDPC, 30 trials, min-sum, max iterations = 20');
xlabel('BER');
ylabel('Error Count');
set(gca, 'xdir', 'reverse');