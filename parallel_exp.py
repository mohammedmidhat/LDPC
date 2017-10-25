from multiprocessing import Pool
import math
import numpy as np
from scipy.stats import bernoulli
import random
import pylab
import sympy
import copy
from scipy import sparse
from copy import deepcopy
from time import clock
import itertools



def parse_MacKay_dense(matrix_file, k, n):
    file_array = np.loadtxt(matrix_file)
    H_matrix = np.array([[0 for j in range(n)] for i in range(n-k)])
    for i in range(n):
        H_matrix[file_array[i,:].astype(int)-1,i] = 1
    return H_matrix


def parse_mine_dense(matrix_file):
    return np.loadtxt(matrix_file)


##def remove_H_redundancy(H):
##    redundant_rows
##    for i in range(len(H)-1):
##        
##    return H


def generate_cw_and_LLRs_run(bit_err_prob, G_t):
    codeword_bases = bernoulli.rvs(0.5, size=G_t.shape[1])
    codeword = np.sum(G_t[:,codeword_bases==1],1)%2
    LLRs = np.array([math.log((1-bit_err_prob)/bit_err_prob) for i in range(G_t.shape[0])])
    noise_word = bernoulli.rvs(bit_err_prob, size=G_t.shape[0])
    LLRs[noise_word+codeword==1] *= -1
    return (LLRs,codeword,noise_word,codeword_bases)


def generate_cw_and_LLRs_AWGN(stdev, G_t):
    codeword_bases = bernoulli.rvs(0.5, size=G_t.shape[1])
    codeword = np.sum(G_t[:,codeword_bases==1],1)%2
    voltages = []
    for i in codeword:
        voltage = np.random.normal(i,stdev,1)[0]
        if(voltage>1):
            voltages.append(2-voltage if i==1 else 1)
        elif(voltage<0):
            voltages.append(-voltage if i==0 else 0)
        else:
            voltages.append(voltage)
    voltages = np.array(voltages)
    
    LLRs = np.array([math.log(math.exp(-i**2/(2*stdev**2))/math.exp(-(i-1)**2/(2*stdev**2))) for i in voltages])
    
    return (LLRs,codeword,voltages)


def LDPC_test_single_trial_AWGN_run(H, G_t, stdev, max_num_iterations):
    err_count = 0
    undet_err_count = 0
    data = generate_cw_and_LLRs_AWGN(stdev, G_t)
    LLRs = data[0]
    cw = data[1]
    decoded_cw = run(H, LLRs, max_num_iterations)
    if(decoded_cw[1]):
        if(not np.all((cw==decoded_cw[0].T[0]))):
            undet_err_count = 1
    else:
        err_count = 1
    return(stdev,err_count,undet_err_count,decoded_cw[2],decoded_cw[1])


def LDPC_test_single_trial_run(H, G_t, bit_err_prob, max_num_iterations):
    err_count = 0
    undet_err_count = 0
    data = generate_cw_and_LLRs_run(bit_err_prob, G_t)
    LLRs = data[0]
    cw = data[1]
    decoded_cw = run(H, LLRs, max_num_iterations)
    if(decoded_cw[1]):
        if(not np.all((cw==decoded_cw[0].T[0]))):
            undet_err_count = 1
    else:
        err_count = 1
    return(bit_err_prob,err_count,undet_err_count,decoded_cw[2],decoded_cw[1])


def LDPC_test_AWGN_run(input_arg):
    H = input_arg[0]
    G_t = input_arg[1]
    stdev = [input_arg[2]]
    max_num_iterations = input_arg[3]
    num_trials = input_arg[4]

    err_count_arr = []
    undet_err_count_arr = []
    iter_arr = []
    hard_success_arr = []
    for i in stdev:
        err_count = 0
        undet_err_count = 0
        iter_count = 0
        hard_success_iter = 0
        for j in range(num_trials):
            result = LDPC_test_single_trial_AWGN_run(H, G_t, i, max_num_iterations)
            err_count += result[1]
            undet_err_count += result[2]
            if(type(result[4])==int):
                iter_count += result[4]
            elif(not result[4]):
                iter_count += max_num_iterations
            else:
                hard_success_iter += 1
        err_count_arr.append(1.0*err_count/num_trials)
        undet_err_count_arr.append(1.0*undet_err_count/num_trials)
        iter_arr.append(iter_count)
        hard_success_arr.append(hard_success_iter)
    return(stdev,err_count_arr,undet_err_count_arr,iter_arr,hard_success_arr)


def LDPC_test_run(input_arg):
    H = input_arg[0]
    G_t = input_arg[1]
    bit_err_prob = [input_arg[2]]
    max_num_iterations = input_arg[3]
    num_trials = input_arg[4]

    err_count_arr = []
    undet_err_count_arr = []
    iter_arr = []
    hard_success_arr = []
    for i in bit_err_prob:
        err_count = 0
        undet_err_count = 0
        iter_count = 0
        hard_success_iter = 0
        for j in range(num_trials):
            result = LDPC_test_single_trial_run(H, G_t, i, max_num_iterations)
            err_count += result[1]
            undet_err_count += result[2]
            if(type(result[4])==int):
                iter_count += result[4]
            elif(not result[4]):
                iter_count += max_num_iterations
            else:
                hard_success_iter += 1
        err_count_arr.append(1.0*err_count/num_trials)
        undet_err_count_arr.append(1.0*undet_err_count/num_trials)
        iter_arr.append(iter_count)
        hard_success_arr.append(hard_success_iter)
    return(bit_err_prob,err_count_arr,undet_err_count_arr,iter_arr,hard_success_arr)


def neighbors_vr(H, vr):
    return H.T[vr].nonzero()[0]

def neighbors_cn(H, cn):
    return H[cn].nonzero()[0]


def LLR_to_cw(LLRs):
    cw = np.array([[0] for i in range(len(LLRs))])
    for i in range(len(LLRs)):
        if(LLRs[i]>0):
            cw[i][0] = 0
        elif(LLRs[i]<0):
            cw[i][0] = 1
    return cw


##def LLR_to_cw(LLRs):
##    cw = np.array([[0] for i in range(len(LLRs))])
##    cw[LLRs<0] = 1
##    return cw

def update_vr_beliefs(H, cn_to_vr_msgs, LLR_init):
    LLRs = np.array([0 for i in range(len(H[0]))])
    for i in range(len(H[0])):
        LLRs[i] = LLR_init[i] + np.sum(cn_to_vr_msgs[neighbors_vr(H,i),i])
    return LLRs


##def update_vr_beliefs(H, cn_to_vr_msgs, LLR_init):
##    LLRs = np.array([0 for i in range(len(H[0]))])
##    for i in range(len(H[0])):
##        LLRs[i] = LLR_init[i] + sum(cn_to_vr_msgs[neighbors_vr(H, i),i])
##    return LLRs


def update_vr_to_cn_msgs(H, cn_to_vr_msgs, vr_to_cn_msgs, LLR_init):
    for i in range(len(H[0])):
        msg = LLR_init[i]
        for j in neighbors_vr(H, i):
            msg += cn_to_vr_msgs[j][i]
        for j in neighbors_vr(H, i):
            vr_to_cn_msgs[i][j] = msg - cn_to_vr_msgs[j][i]


##def update_vr_to_cn_msgs(H, cn_to_vr_msgs, vr_to_cn_msgs, LLRs):
##    for i in range(len(H[0])):
##        vr_to_cn_msgs[i,neighbors_vr(H, i)] = LLRs[i] - cn_to_vr_msgs[neighbors_vr(H, i),i]


def update_cn_to_vr_msgs(H, vr_to_cn_msgs, cn_to_vr_msgs, LLRs):
    for j in range(len(H)):
        #print('ck11')
        for i in neighbors_cn(H, j):
            #print('ck12')
            product = 1.0
            for k in neighbors_cn(H, j):
                #print('ck13')
                if(k!=i):
                    product*=math.tanh(0.5*vr_to_cn_msgs[k][j])
                    #print('ck14')
            #print('ck15')
            if(abs(product)<0.9999999999999999):
                cn_to_vr_msgs[j][i] = 2*math.atanh(product)
            else:
                cn_to_vr_msgs[j][i] = 2*math.atanh(0.9999999999999999)
            #print('ck16')


##def update_cn_to_vr_msgs(H, vr_to_cn_msgs, cn_to_vr_msgs, LLRs):
##    for j in range(len(H)):
##        #print('ck11')
##        for i in neighbors_cn(H, j):
##            #print('ck12')
##            msg_list = []
##            for k in neighbors_cn(H, j):
##                #print('ck13')
##                if(k!=i):
##                    msg_list.append(vr_to_cn_msgs[k][j])
##                    #print('ck14')
##            #print('ck15')
##            msg_list = np.array(msg_list)
##            cn_to_vr_msgs[j][i] = np.prod(np.sign(msg_list))*np.min(abs(msg_list))
            


##def update_cn_to_vr_msgs(H, vr_to_cn_msgs, cn_to_vr_msgs):
##    for j in range(len(H)):
##        #print('ck11')
##        for i in neighbors_cn(H, j):
##            #print('ck12')
##            neighbors_in_prod = neighbors_cn(H, j)[neighbors_cn(H, j)!=i]
##            product=np.prod(np.tanh(0.5*vr_to_cn_msgs[neighbors_in_prod,j]))
##            #print('ck15')
##            if(abs(product)<0.9999999999999999):
##                cn_to_vr_msgs[j][i] = 2*math.atanh(product)
##            else:
##                cn_to_vr_msgs[j][i] = 2*math.atanh(0.9999999999999999)
##            #print('ck16')


def vr_to_cn_msgs_init(LLR_init, H):
    result = []
    for i in range(len(H[0])):
        row_in_result = [0 for j in range(len(H))]
        for j in neighbors_vr(H, i):
            row_in_result[j] = LLR_init[i]
        result.append(row_in_result)
    return np.array(result)


##def vr_to_cn_msgs_init(LLR_init, H):
##    result = np.array([[0 for j in range(len(H))] for i in range(len(H[0]))])
##    for i in range(len(H[0])):
##        result[i,neighbors_vr(H, i)] = LLR_init[i]
##    return result


def run(H, LLR_init, max_iterations):
    t1=clock()
    if(np.array_equal(np.dot(H,LLR_to_cw(LLR_init))%2, np.array([[0] for j in range(len(H))]))):
        #print('ck7')
        return (LLR_to_cw(LLR_init)%2, True, clock()-t1)
    #print('ck1')
    vr_to_cn_msgs = vr_to_cn_msgs_init(LLR_init, H)
    #print('ck2')
    cn_to_vr_msgs = np.array([[0 for j in range(len(H[0]))] for i in range(len(H))])
    #print('ck3')
    for i in range(max_iterations):
        #print('ck4')
        update_cn_to_vr_msgs(H, vr_to_cn_msgs, cn_to_vr_msgs, LLR_init)
        #print('ck5')
        LLRs = update_vr_beliefs(H, cn_to_vr_msgs, LLR_init)
        #print('ck6')
        if(np.array_equal(np.dot(H,LLR_to_cw(LLRs))%2, np.array([[0] for j in range(len(H))]))):
            #print('ck7')
            return (LLR_to_cw(LLRs)%2, i+1, clock()-t1)
        update_vr_to_cn_msgs(H, cn_to_vr_msgs, vr_to_cn_msgs, LLR_init)
        #print('ck8')
    return (LLR_to_cw(LLRs)%2, False, clock()-t1)


def graph_SFR_BSC(bit_err_prob,err_count,undet_err_count):
    pylab.plot(bit_err_prob,err_count, '-o', label = 'Error Count')
    pylab.plot(bit_err_prob,undet_err_count,  '-o', label = 'Undet Error Count')
    pylab.gca().invert_xaxis()
    pylab.grid(True)
    pylab.xlabel('BER')
    pylab.legend(loc=1)
    pylab.show()


def graph_SFR_BSC_log_scale(bit_err_prob,err_count,undet_err_count):
    pylab.semilogy(bit_err_prob,err_count, '-o', label = 'Error Count')
    pylab.semilogy(bit_err_prob,undet_err_count,  '-o', label = 'Undet Error Count')
    pylab.gca().invert_xaxis()
    pylab.grid(True)
    pylab.xlabel('BER')
    pylab.legend(loc=1)
    pylab.show()


def graph_SFR_AWGN(stdev,err_count,undet_err_count):
    pylab.plot(stdev,err_count, '-o', label = 'Error Count')
    pylab.plot(stdev,undet_err_count,  '-o', label = 'Undet Error Count')
    pylab.gca().invert_xaxis()
    pylab.grid(True)
    pylab.xlabel('stdev')
    pylab.legend(loc=1)
    pylab.show()


def graph_SFR_AWGN_SNR(stdev,err_count,undet_err_count):
    SNR_in_dB = 10*np.log10(1/np.square(stdev))
    pylab.plot(SNR_in_dB,err_count, '-o', label = 'Error Count')
    pylab.plot(SNR_in_dB,undet_err_count,  '-o', label = 'Undet Error Count')
    pylab.grid(True)
    pylab.xlabel('SNR (dB)')
    pylab.legend(loc=1)
    pylab.show()


def graph_SFR_AWGN_log_scale(stdev,err_count,undet_err_count):
    pylab.semilogy(stdev,err_count, '-o', label = 'Error Count')
    pylab.semilogy(stdev,undet_err_count,  '-o', label = 'Undet Error Count')
    pylab.gca().invert_xaxis()
    pylab.grid(True)
    pylab.xlabel('stdev')
    pylab.legend(loc=1)
    pylab.show()


def graph_SFR_AWGN_log_scale_SNR(stdev,err_count,undet_err_count):
    SNR_in_dB = 10*np.log10(1/np.square(stdev))
    pylab.semilogy(SNR_in_dB,err_count, '-o', label = 'Error Count')
    pylab.semilogy(SNR_in_dB,undet_err_count,  '-o', label = 'Undet Error Count')
    pylab.grid(True)
    pylab.xlabel('SNR (dB)')
    pylab.legend(loc=1)
    pylab.show()

##H = parse_mine_dense('H_proto_based_1.txt')
##G_t = np.loadtxt('G_t_proto_based_1.txt')
##H = parse_MacKay_dense('204.33.484 (N=204,K=102,M=102,R=0.5).txt',102,204)
##G_t = np.loadtxt('204.33.484 (N=204,K=102,M=102,R=0.5)G_t.txt')
##t = clock()
##result = LDPC_test_run([H, G_t, 0.05, 10, 1])
##result = LDPC_test_AWGN_run([H, G_t, 0.05, 10, 1])
##t_post = clock()-t

if __name__ == '__main__':
    H = parse_mine_dense('H_proto_based_0.txt')
    G_t = np.loadtxt('G_t_proto_based_0.txt')
    #H = parse_MacKay_dense('204.33.484 (N=204,K=102,M=102,R=0.5).txt',102,204)
    #G_t = np.loadtxt('204.33.484 (N=204,K=102,M=102,R=0.5)G_t.txt')
    p = Pool(4)
    t = clock()
    result = p.map(LDPC_test_AWGN_run, [(H,G_t,0.01*i,10,10000) for i in range(30,3,-1)])
    #result = p.map(LDPC_test_run, [(H,G_t,0.01*i,10,1000) for i in range(10,1,-1)])
    t_post = clock()-t
    graph_SFR_AWGN_log_scale_SNR([i[0][0] for i in result],[i[1][0] for i in result],[i[2][0] for i in result])
    #graph_SFR_BSC_log_scale([i[0][0] for i in result],[i[1][0] for i in result],[i[2][0] for i in result])
