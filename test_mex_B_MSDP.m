conf_mat = 0.95*eye(16);
conf_mat(2,1) = 0.05;
conf_mat(12,13) = 0.05;
for i = 2:12
    conf_mat(i-1,i) = 0.025;
    conf_mat(i+1,i) = 0.025;
end

bit_err_prob = [0.025:-0.005:0.001];
err_count = zeros(1,size(bit_err_prob,2));
undet_err = zeros(1,size(bit_err_prob,2));
num_trials = 4;
max_num_iter = 30;
num_reads = 1;

result = test_code_B_MSDP(max_num_iter,num_trials, num_reads, 1, conf_mat');

% tic;
% parfor (j = 1:size(bit_err_prob,2),4)
%     result = test_code_B(max_num_iter,num_trials,bit_err_prob(1,j));
%     %result = test_code_min_sum_B(max_num_iter,num_trials,bit_err_prob(1,j));
%     err_count(1,j) = result(1)/num_trials;
%     undet_err(1,j) = result(2)/num_trials;
% end
% t_post = toc;
% 
% loglog(bit_err_prob, err_count,'--o');
% hold on;
% loglog(bit_err_prob, undet_err,'--o');
% title('(64800,48600) LDPC, INT = 8 bits, FRAC = 14 bits, max iterations = 30, dec_alg = SPA, DVB');
% xlabel('Bit Error Rate (BER)');
% ylabel('Frame Error Rate (FER)');
% legend('error','undetected error');
% set(gca, 'xdir', 'reverse');
