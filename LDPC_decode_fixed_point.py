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
from deModel import DeFixedInt as fix

int_width = 16
frac_width = 4

def get_fixed(x):
    return fix(int_width,frac_width,x)._getFloatValue()

def parse_MacKay_dense_fix(matrix_file, k, n):
    file_array = np.loadtxt(matrix_file)
    H_matrix = np.array([[0 for j in range(n)] for i in range(n-k)])
    for i in range(n):
        H_matrix[file_array[i,:].astype(int)-1,i] = 1
    return H_matrix

def generate_cw_and_LLRs_run_fix(bit_err_prob, G_t):
    codeword_bases = bernoulli.rvs(0.5, size=G_t.shape[1])
    codeword = np.sum(G_t[:,codeword_bases==1],1)%2
    LLRs = np.array([get_fixed(math.log((1-bit_err_prob)/bit_err_prob)) for i in range(G_t.shape[0])])
    noise_word = bernoulli.rvs(bit_err_prob, size=G_t.shape[0])
    LLRs[noise_word+codeword==1] *= -1
    return (LLRs,codeword,noise_word,codeword_bases)


def LDPC_test_single_trial_run_fix(H, G_t, bit_err_prob, max_num_iterations):
    err_count = 0
    undet_err_count = 0
    data = generate_cw_and_LLRs_run_fix(bit_err_prob, G_t)
    LLRs = data[0]
    cw = data[1]
    decoded_cw = run_fix(H, LLRs, max_num_iterations)
    if(decoded_cw[1]):
        if(not np.all((cw==decoded_cw[0].T[0]))):
            undet_err_count = 1
    else:
        err_count = 1
    return(bit_err_prob,err_count,undet_err_count,decoded_cw[2])


def LDPC_test_run_fix(input_arg):
    H = input_arg[0]
    G_t = input_arg[1]
    bit_err_prob = [input_arg[2]]
    max_num_iterations = input_arg[3]
    num_trials = input_arg[4]
    err_count_arr = []
    undet_err_count_arr = []
    for i in bit_err_prob:
        err_count = 0
        undet_err_count = 0
        for j in range(num_trials):
            result = LDPC_test_single_trial_run_fix(H, G_t, i, max_num_iterations)
            err_count += result[1]
            undet_err_count += result[2]
        err_count_arr.append(1.0*err_count/num_trials)
        undet_err_count_arr.append(1.0*undet_err_count/num_trials)
    return(bit_err_prob,err_count_arr,undet_err_count_arr)


def neighbors_vr_fix(H, vr):
    return H.T[vr].nonzero()[0]

def neighbors_cn_fix(H, cn):
    return H[cn].nonzero()[0]


def LLR_to_cw_fix(LLRs):
    cw_string = ''
    cw = np.array([[0] for i in range(len(LLRs))])
    for i in range(len(LLRs)):
        if(LLRs[i]>0):
            cw[i][0] = 0
            cw_string += '0'
        elif(LLRs[i]<0):
            cw[i][0] = 1
            cw_string += '1'
    return (cw,cw_string)

def update_vr_beliefs_fix(H, cn_to_vr_msgs, LLR_init):
    LLRs = []
    for i in range(len(H[0])):
        belief = LLR_init[i]
        for j in neighbors_vr_fix(H, i):
            belief += cn_to_vr_msgs[j][i]
        LLRs.append(belief)
    return LLRs

def update_vr_to_cn_msgs_fix(H, cn_to_vr_msgs, vr_to_cn_msgs, LLR_init):
    for i in range(len(H[0])):
        msg = LLR_init[i]
        for j in neighbors_vr_fix(H, i):
            msg += cn_to_vr_msgs[j][i]
        for j in neighbors_vr_fix(H, i):
            vr_to_cn_msgs[i][j] = msg - cn_to_vr_msgs[j][i]


def update_cn_to_vr_msgs_min_fix(H, vr_to_cn_msgs, cn_to_vr_msgs, LLRs):
    for j in range(len(H)):
        #print('ck11')
        for i in neighbors_cn_fix(H, j):
            #print('ck12')
            msg_list = []
            for k in neighbors_cn_fix(H, j):
                #print('ck13')
                if(k!=i):
                    msg_list.append(vr_to_cn_msgs[k][j])
                    #print('ck14')
            #print('ck15')
            msg_list = np.array(msg_list)
            cn_to_vr_msgs[j][i] = np.prod(np.sign(msg_list))*np.min(abs(msg_list))


##def update_cn_to_vr_msgs_fix(H, vr_to_cn_msgs, cn_to_vr_msgs, LLRs):
##    for j in range(len(H)):
##        #print('ck11')
##        for i in neighbors_cn_fix(H, j):
##            #print('ck12')
##            product = 1.0
##            for k in neighbors_cn_fix(H, j):
##                #print('ck13')
##                if(k!=i):
##                    product*=get_fixed(math.tanh(0.5*vr_to_cn_msgs[k][j]))
##                    #print('ck14')
##            #print('ck15')
##            if(abs(product)<0.9999999999999999):
##                cn_to_vr_msgs[j][i] = get_fixed(2*get_fixed(math.atanh(product)))
##            else:
##                cn_to_vr_msgs[j][i] = get_fixed(2*get_fixed(math.atanh(0.9999999999999999)))
##            #print('ck16')


def vr_to_cn_msgs_init_fix(LLR_init, H):
    result = []
    for i in range(len(H[0])):
        row_in_result = [0 for j in range(len(H))]
        for j in neighbors_vr_fix(H, i):
            row_in_result[j] = LLR_init[i]
        result.append(row_in_result)
    return result
    

def run_fix(H, LLR_init, max_iterations):
    #print('ck1')
    vr_to_cn_msgs = vr_to_cn_msgs_init_fix(LLR_init, H)
    #print('ck2')
    cn_to_vr_msgs = [[0 for j in range(len(H[0]))] for i in range(len(H))]
    #print('ck3')
    t1=clock()
    for i in range(max_iterations):
        #print vr_to_cn_msgs
        #print('ck4')
        update_cn_to_vr_msgs_min_fix(H, vr_to_cn_msgs, cn_to_vr_msgs, LLR_init)
        #print('ck5')
        print cn_to_vr_msgs
        LLRs = update_vr_beliefs_fix(H, cn_to_vr_msgs, LLR_init)
        #print LLRs
        #print('ck6')
        if(np.array_equal(np.dot(H,LLR_to_cw_fix(LLRs)[0])%2, np.array([[0] for j in range(len(H))]))):
            #print('ck7')
            return (LLR_to_cw_fix(LLRs)[0]%2, i+1, clock()-t1)
        update_vr_to_cn_msgs_fix(H, cn_to_vr_msgs, vr_to_cn_msgs, LLR_init)
        #print('ck8')
    return (LLR_to_cw_fix(LLRs)[0]%2, False, clock()-t1)


def graph_SFR(bit_err_prob,err_count,undet_err_count):
    pylab.plot(bit_err_prob,err_count, '-o', label = 'Error Count')
    pylab.plot(bit_err_prob,undet_err_count,  '-o', label = 'Undet Error Count')
    pylab.gca().invert_xaxis()
    pylab.grid(True)
    pylab.xlabel('BER')
    pylab.legend(loc=1)
    pylab.show()

    
##H = parse_MacKay_dense_fix('204.33.484 (N=204,K=102,M=102,R=0.5).txt',102,204)
##G_t = np.loadtxt('204.33.484 (N=204,K=102,M=102,R=0.5)G_t.txt')
##t = clock()
##result = LDPC_test_run_fix([H, G_t, 0.01, 3, 3])
##t_post = clock()-t


## For FPGA testing
##H=np.array([[1,1,1,1,1,1],[1,1,1,1,1,1],[1,1,1,1,1,1]])
##LLR_init = np.array([get_fixed(4.0) for i in range(6)])
##LLR_init[4] = get_fixed(-4.0)
##x=run_fix(H, LLR_init, 6)

##if __name__ == '__main__':
##    ##H = parse_MacKay_dense_fix('204.33.484 (N=204,K=102,M=102,R=0.5).txt',102,204)
##    ##G_t = np.loadtxt('204.33.484 (N=204,K=102,M=102,R=0.5)G_t.txt')
##    p = Pool(4)
##    t = clock()
##    result = p.map(LDPC_test_run_fix, [(H,G_t,0.01*i,10,100) for i in range(5,1,-1)])
##    t_post = clock()-t

H = np.array([[1,1,0,1,0,0],[0,1,1,0,1,0],[1,0,1,0,0,1]])
G_T = np.array([[1,0,1],[1,1,0],[1,0,0],[0,1,1],[0,1,0],[0,0,1]])

