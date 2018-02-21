% This script invokes binary LDPC decoding simulations.
% A confusion matrix characterizes the noise of the channel
% The decoding is invoked through the test_code_B_MSDP() mex call.

% Author: Mohammed Al Ai Baky
% Created: 1/20/2018


conf_mat = 0.948*eye(16);
conf_mat(2,1) = 0.052;
conf_mat(12,13) = 0.052;
for i = 2:12
    conf_mat(i-1,i) = 0.026;
    conf_mat(i+1,i) = 0.026;
end

conf_mat = [ 0.9605,  0.0205,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    
        0.0395,  0.9606,  0.0189,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    
        0.    ,  0.0189,  0.9621,  0.0189,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    
        0.    ,  0.    ,  0.0189,  0.9621,  0.0189,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    
        0.    ,  0.    ,  0.    ,  0.0189,  0.9621,  0.0189,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    
        0.    ,  0.    ,  0.    ,  0.    ,  0.0189,  0.9621,  0.0189,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    
        0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.0189,  0.9621,  0.0189,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    
        0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.0189,  0.9621,  0.0189,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    
        0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.0189,  0.9621,  0.0189,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    
        0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.0189,  0.9621,  0.0189,  0.    ,  0.    ,  0.    ,  0.    ,  0.    
        0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.0189,  0.9621,  0.0189,  0.    ,  0.    ,  0.    ,  0.    
        0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.0189,  0.9621,  0.0189,  0.    ,  0.    ,  0.    
        0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.0189,  0.9621,  0.0189,  0.    ,  0.    
        0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.0189,  0.9621,  0.0189,  0.    
        0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.0189,  0.9621,  0.0189
        0.    ,  0.    , -0.    , -0.    , -0.    , -0.    , -0.    , -0.    , -0.    , -0.    , -0.    , -0.    , -0.    ,  0.    ,  0.0189,  0.9811];

num_trials = 1000;
max_num_iter = 20;
num_reads = 1;
decode_mode = 0;

result = test_code_B_MSDP(max_num_iter,num_trials, num_reads, decode_mode, conf_mat');

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
