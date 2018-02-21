% This script invokes binary LDPC decoding simulations 
% The decoding is invoked through the test_code_B_float() mex call.

% Author: Mohammed Al Ai Baky
% Created: 1/20/2018

cw_noise_gen = 0;

bit_err_prob = [0.1:-0.01:0.04];
total_err_count = cell(1,size(bit_err_prob,2));
total_undet_err = cell(1,size(bit_err_prob,2));
num_trials = 10000;
serial = 100;
parallel = num_trials/serial;
max_num_iter = 30;

tic;

if cw_noise_gen
    for j = 1:size(bit_err_prob,2)
        err_count = zeros(1,parallel);
        undet_err = zeros(1,parallel);
        
        parfor (k = 1:parallel)
            result = test_code_B_float('Frolov_1024_0.5.txt', max_num_iter,serial,bit_err_prob(1,j));
            %result = test_code_min_sum_B('204.33.484.txt', max_num_iter,num_trials,bit_err_prob(1,j));
            err_count(1,k) = result(1)/num_trials;
            undet_err(1,k) = result(2)/num_trials;
        end
        
        total_err_count{j} = [err_count];
        total_undet_err{j} = [undet_err];
        
        save('full_env.mat');
    end
else
    [n, m, cmax, rmax, col_weight, row_weight, row_list_col, col_list_row, col_list_N] = initdec('Frolov_1024_0.5.txt');
    for j = 1:size(bit_err_prob,2)
        err_count = zeros(1,parallel);
        undet_err = zeros(1,parallel);
        
        parfor (k = 1:parallel)
            for i = 1:serial
                % genrate n-vector of random data
                x = randi([0 1],1,n);
                % Compute syndrome
                s = enc(x, row_weight, m, row_list_col);
                % Add noise
                noise = rand(1, n)<bit_err_prob(1,j);
                y = mod(x+noise,2);
                % Compute LLR's
                q0 = zeros(1,n);
                q0(y == 0) = log((1-bit_err_prob(1,j))/bit_err_prob(1,j));
                q0(y == 1) = log(bit_err_prob(1,j)/(1-bit_err_prob(1,j)));
                % Decode
                result = test_code_B_float('Frolov_1024_0.5.txt', max_num_iter, n, m, x, s, y, q0);
                %result = test_code_min_sum_B('204.33.484.txt', max_num_iter, n, m, x, s, y, q0);
                % Collect Results
                err_count(1,k) = result(1)/num_trials;
                undet_err(1,k) = result(2)/num_trials;
            end
        end
        total_err_count{j} = [err_count];
        total_undet_err{j} = [undet_err];
        
        save('full_env.mat');
    end
end

t_post = toc

err_rate = zeros(1,size(bit_err_prob,2));
undet_err_rate = zeros(1,size(bit_err_prob,2));
for i = 1:size(bit_err_prob,2)
    err_rate(1,i) = sum(total_err_count{i});
    undet_err_rate(1,i) = sum(total_undet_err{i});
end

loglog(bit_err_prob, err_rate,'--o');
%hold on;
loglog(bit_err_prob, undet_err_rate,'--o');
title('(2048,1024) LDPC, max iterations = 30, dec alg = SPA, Frolov');
xlabel('Bit Error Rate (BER)');
ylabel('Frame Error Rate (FER)');
legend('error','undetected error');
set(gca, 'xdir', 'reverse');
