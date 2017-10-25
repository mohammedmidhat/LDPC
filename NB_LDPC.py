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
from multiprocessing import Pool
import itertools


#np.set_printoptions(threshold=np.nan)
inverse_q4 = {1:1,2:3,3:2}
#inverse_q8 = {1:1,2:,3:,4:,5:,6:,7:}

def multip_q4(M,x):
    new_mat = M*x
    new_mat[new_mat==9] = 2
    new_mat[new_mat==6] = 1
    new_mat[new_mat==4] = 3
    return np.sum(new_mat,1)%4


def parse_MacKay_dense(matrix_file, k, n):
    file_array = np.loadtxt(matrix_file)
    H_matrix = np.array([[0 for j in range(n)] for i in range(n-k)])
    for i in range(n):
        H_matrix[file_array[i,:].astype(int)-1,i] = 1
    return H_matrix


def parse_MacKay_NB(matrix_file, k, n):
    file_array = np.loadtxt(matrix_file)
    H_matrix = np.array([[0 for j in range(n)] for i in range(n-k)])
    for i in range(n):
        H_matrix[file_array[i,[0,2,4,6]].astype(int)-1,i] = file_array[i,[1,3,5,7]]
    return H_matrix


def generate_cw_and_LLRs_AWGN_run_NB(stdev, G_t, q):
    user_data = np.random.random_integers(0, q-1, size=(1,G_t.shape[1]))
    codeword = multip_q4(G_t,user_data)

    voltages = []

    for i in codeword:
        voltage = np.random.normal(1.0*i/(q-1),stdev,1)[0]
        if(voltage>1):
            voltages.append(2-voltage if i==q-1 else 1)
        elif(voltage<0):
            voltages.append(-voltage if i==0 else 0)
        else:
            voltages.append(voltage)
    voltages = np.array(voltages)
    LLRs = np.array([[math.log(math.exp(-(i-1.0*j/(q-1))**2/(2*stdev**2))/math.exp(-i**2/(2*stdev**2))) for j in range(1,q)] for i in voltages])
    return (LLRs, codeword)



def LDPC_test_single_trial_AWGN_run_NB(H, G_t, stdev, max_num_iterations):
    undet_err_cws = []
    cws_all = []
    err_count = 0
    undet_err_count = 0
    data = generate_cw_and_LLRs_AWGN_run_NB(stdev, G_t)
    LLRs = data[0]
    cw = data[1]
    decoded_cw = run_NB(H, LLRs, max_num_iterations)
    if(decoded_cw[1]):
        if(not np.all((cw==decoded_cw[0].T[0]))):
            undet_err_count = 1
            undet_err_cws.append((cw,decoded_cw[0]))
    else:
        err_count = 1
    cws_all.append((cw,decoded_cw[0]))
    return(stdev,err_count,undet_err_count,decoded_cw[2],undet_err_cws,cws_all,LLRs)


def LDPC_test_AWGN_run_NB(H, G_t, bit_err_prob, max_num_iterations, num_trials):
    LLRs = []
    undet_err_cws = []
    cws_all = []
    err_count_arr = []
    undet_err_count_arr = []
    for i in bit_err_prob:
        err_count = 0
        undet_err_count = 0
        for j in range(num_trials):
            result = LDPC_test_single_trial_run(H, G_t, i, max_num_iterations)
            err_count += result[1]
            undet_err_count += result[2]
            undet_err_cws.append(result[4])
            cws_all.append(result[5])
            LLRs.append(result[6])
        err_count_arr.append(err_count)
        undet_err_count_arr.append(undet_err_count)
    return(bit_err_prob,err_count_arr,undet_err_count_arr,undet_err_cws,cws_all,LLRs)

#########################################################################
#########################################################################
#########################################################################
def permuter_q4(vec,perm):
    if(perm==2):
        return vec[[1,2,0]]
    elif(perm==3):
        return vec[[2,0,1]]
    else:
        return vec


def perm_q4(perm):
    if(perm==2):
        return [1,2,0]
    elif(perm==3):
        return [2,0,1]
    else:
        return [0,1,2]
    

def build_P_nb(H_nb,q):
    P = np.array([[[0 for k in range(q-1)] for j in range(len(H_nb[0]))] for i in range(len(H_nb))])
    for i in range(1,q):
        P[H_nb==i] = perm_q4(i)
    return P


def neighbors_vr(H, vr):
    return H.T[vr].nonzero()[0]


def neighbors_cn(H, cn):
    return H[cn].nonzero()[0]


def LLR_to_cw(LLRs):
    cw = np.array([[0] for i in range(len(LLRs))])
    for i in range(len(LLRs)):
        if(np.max(LLRs[i]) > 0):
            cw[i][0] = np.argmax(LLRs[i])+1
        else:
            cw[i][0] = 0
    return cw


def update_p_to_vr_msgs(H_nb, P_nb, p_to_vr_msgs, cn_to_p_msgs):
    for i in range(len(H_nb)):
        p_to_vr_msgs[i,neighbors_cn(H_nb,i),:] = cn_to_p_msgs[i,neighbors_cn(H_nb,i), P_nb[i,neighbors_cn(H_nb,i)]]
    

def update_vr_beliefs(H, p_to_vr_msgs, LLR_init):
    LLRs = []
    for i in range(len(H[0])):
        belief = LLR_init[i]
        for j in neighbors_vr(H, i):
            belief += cn_to_vr_msgs[j][i]
        LLRs.append(belief)
    return LLRs


def update_vr_to_p_msgs(H, p_to_vr_msgs, vr_to_p_msgs, LLR_init):
    for i in range(len(H[0])):
        msg = LLR_init[i]
        for j in neighbors_vr(H, i):
            msg += p_to_vr_msgs[j][i]
        for j in neighbors_vr(H, i):
            p_to_cn_msgs[i][j] = msg - p_to_vr_msgs[j][i]


def update_cn_to_vr_msgs(H, vr_to_cn_msgs, cn_to_vr_msgs, LLRs):
    for j in range(len(H)):
        #print('ck11')
        for i in neighbors_cn(H, j):
            #print('ck12')
            product = 1.0
            #print('neighbor ' + str(i) + 'of cn '+ str(j))
            for k in neighbors_cn(H, j):
                #print('ck13')
                if(k!=i):
                    product*=math.tanh(0.5*vr_to_cn_msgs[k][j])
                    #print('ck14')
                #print('product = ' + str(product))
            #print('ck15')
            if(abs(product)<0.9999999999999999):
                #print('less')
                cn_to_vr_msgs[j][i] = 2*math.atanh(product)
            else:
                #print('more')
                cn_to_vr_msgs[j][i] = 2*math.atanh(0.9999999999999999)
            #print('ck16')


def vr_to_cn_msgs_init(LLR_init, H):
    result = []
    for i in range(len(H[0])):
        row_in_result = [0 for j in range(len(H))]
        for j in neighbors_vr(H, i):
            row_in_result[j] = LLR_init[i]
        result.append(row_in_result)
    return result
    

def run(H, LLR_init, max_iterations):
    #print('ck1')
    vr_to_cn_msgs = vr_to_cn_msgs_init(LLR_init, H)
    #print('ck2')
    cn_to_vr_msgs = [[0 for j in range(len(H[0]))] for i in range(len(H))]
    #print('ck3')
    t1=clock()
    for i in range(max_iterations):
        #print('ck4')
        update_cn_to_vr_msgs(H, vr_to_cn_msgs, cn_to_vr_msgs, LLR_init)
        #print('ck5')
        LLRs = update_vr_beliefs(H, cn_to_vr_msgs, LLR_init)
        #print('ck6')
        if(np.array_equal(np.dot(H,LLR_to_cw(LLRs)[0])%2, np.array([[0] for j in range(len(H))]))):
            #print('ck7')
            return (LLR_to_cw(LLRs)[0]%2, i+1, clock()-t1)
        update_vr_to_cn_msgs(H, cn_to_vr_msgs, vr_to_cn_msgs, LLR_init)
        #print('ck8')
    return (LLR_to_cw(LLRs)[0]%2, False, clock()-t1)


##H = parse_MacKay_dense('204.33.484 (N=204,K=102,M=102,R=0.5).txt',102,204)
##G_t = np.loadtxt('204.33.484 (N=204,K=102,M=102,R=0.5)G_t.txt')


##if __name__ == '__main__':
##    p = Pool(4)
##    print(p.map(randperm, [10,17, 7, 9]))
    
