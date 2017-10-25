def write_matlab_script(filename,BER,):
    


processors_per_comp = 32
bit_err_prob = [i/10.0 for i in range(1,11)]
num_serial_trials = 1000
num_trials = 1000000

num_parallel_trials = 1.0*num_trials/num_serial_trials

for i in bit_err_prob:
    for j in range(num_parallel_trials):
        
