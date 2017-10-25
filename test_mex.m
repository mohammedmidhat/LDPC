bit_err_prob = [0.2:-0.01:0.05];
err_count = zeros(1,size(bit_err_prob,2));
undet_err = zeros(1,size(bit_err_prob,2));
num_trials = 1000;
max_num_iter = 30;

tic;
parfor (j = 1:size(bit_err_prob,2),4)
    result = test_code_B(max_num_iter,num_trials,bit_err_prob(1,j));
    %result = test_code_min_sum_B(max_num_iter,num_trials,bit_err_prob(1,j));
    display('no');
    err_count(1,j) = result(1)/num_trials;
    undet_err(1,j) = result(2)/num_trials;
end
t_post = toc;

loglog(bit_err_prob, err_count,'--o');
hold on;
loglog(bit_err_prob, undet_err,'--o');
title('(64800,48600) LDPC, INT = 8 bits, FRAC = 14 bits, max iterations = 30, dec_alg = SPA, DVB');
xlabel('Bit Error Rate (BER)');
ylabel('Frame Error Rate (FER)');
legend('error','undetected error');
set(gca, 'xdir', 'reverse');
