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


np.set_printoptions(threshold=np.nan)


# John Fan's book example
#H = np.array([[1,1,0,1,0,0],[0,1,1,0,1,0],[1,0,1,0,0,1]])


def parse_MacKay_dense(matrix_file, k, n):
    file_array = np.loadtxt(matrix_file)
    H_matrix = np.array([[0 for j in range(n)] for i in range(n-k)])
    for i in range(n):
        H_matrix[file_array[i,:].astype(int)-1,i] = 1
    return H_matrix


def parse_MacKay_sparse(matrix_file, k, n):
    file_array = np.loadtxt(matrix_file)
    H_matrix = sparse.lil_matrix([[0 for j in range(n)] for i in range(n-k)])
    for i in range(n):
        H_matrix[file_array[i,:]-1,i] = 1
    return H_matrix


def decode(H, LLRs, max_num_iterations):
    nCheck = H.shape[0]
    nBit = H.shape[1]
    indexCheck = H.T.nonzero()[1]
    indexBit = H.T.nonzero()[0]
    indexBit_H = H.nonzero()[1]
    indexCheck_H = H.nonzero()[0]
    cw = sparse.lil_matrix([0 for i in range(nBit)])
##    checks = [0 for i in range(nCheck)]
##    for i in range(nCheck):
##        checks[i] = indexBit[indexCheck==i]
    messageBit = LLRs
    msgb2ch = sparse.lil_matrix([[0 for j in range(nBit)] for i in range(nCheck)], dtype = 'float64')
    msgch2b = sparse.lil_matrix([[0 for j in range(nBit)] for i in range(nCheck)], dtype = 'float64')
    t1 = clock()
    for i in range(max_num_iterations):
        for j in range(nCheck):
            msgb2ch[j,indexBit[indexCheck==j]] = messageBit[0,indexBit[indexCheck==j]] - msgch2b[j,indexBit[indexCheck==j]]
            # indexBit[indexCheck==j] gives the neighbors of j
            for k in range(len(indexBit[indexCheck==j])):
                prod = np.prod(np.tanh((msgb2ch[j,np.delete(indexBit[indexCheck==j],k)]/2.0).toarray()))
                if(abs(prod)>0.99):
                    prod = 0.9999999999999999
                msgch2b[j,indexBit[indexCheck==j][k]] = 2*math.atanh(prod)
        messageBit = sum(msgch2b) + LLRs
        cw[messageBit > 0] = 0
        cw[messageBit < 0] = 1
        if(np.count_nonzero((H*cw.T).toarray()%2)==0):
            return (i+1, cw, clock()-t1)          ## return (num_iterations, cw)
    return (False,0,clock()-t1)


# returns   LLRs: sparse.lil_matrix
#           codeword: sparse.lil_matrix
#           noise_word: sparse.lil_matrix
def generate_cw_and_LLRs(bit_err_prob, G_t):
    codeword_bases = bernoulli.rvs(0.5, size=G_t.shape[1])
    codeword = np.sum(G_t[:,codeword_bases==1],1)%2
    LLRs = sparse.lil_matrix([math.log((1-bit_err_prob)/bit_err_prob) for i in range(G_t.shape[0])])
    noise_word = bernoulli.rvs(bit_err_prob, size=G_t.shape[0])
    LLRs[0, noise_word+codeword==1] *= -1
    return (LLRs,sparse.lil_matrix(codeword),sparse.lil_matrix(noise_word))


def LDPC_test_single_trial(H, G_t, bit_err_prob, max_num_iterations):
    err_count = 0
    undet_err_count = 0
    data = generate_cw_and_LLRs(bit_err_prob, G_t)
    LLRs = data[0]
    cw = data[1]
    decoded_cw = decode(H, LLRs, max_num_iterations)
    if(decoded_cw[0]):
        if(not np.all((cw==decoded_cw[1]).toarray())):
            undet_err_count = 1
    else:
        err_count = 1
    return(bit_err_prob,err_count,undet_err_count,decoded_cw[2])


def LDPC_test(H, G_t, bit_err_prob, max_num_iterations, num_trials):
    err_count_arr = []
    undet_err_count_arr = []
    for i in bit_err_prob:
        err_count = 0
        undet_err_count = 0
        for j in range(num_trials):
            result = LDPC_test_single_trial(H, G_t, i, max_num_iterations)
            err_count += result[1]
            undet_err_count += result[2]
        err_count_arr.append(err_count)
        undet_err_count_arr.append(undet_err_count)
    return(bit_err_prob,err_count_arr,undet_err_count_arr)


def generate_cw_and_LLRs_run(bit_err_prob, G_t):
    codeword_bases = bernoulli.rvs(0.5, size=G_t.shape[1])
    codeword = np.sum(G_t[:,codeword_bases==1],1)%2
    LLRs = np.array([math.log((1-bit_err_prob)/bit_err_prob) for i in range(G_t.shape[0])])
    noise_word = bernoulli.rvs(bit_err_prob, size=G_t.shape[0])
    LLRs[noise_word+codeword==1] *= -1
    return (LLRs,codeword,noise_word,codeword_bases)


def LDPC_test_single_trial_run(H, G_t, bit_err_prob, max_num_iterations):
    undet_err_cws = []
    cws_all = []
    err_count = 0
    undet_err_count = 0
    data = generate_cw_and_LLRs_run(bit_err_prob, G_t)
    LLRs = data[0]
    cw = data[1]
    decoded_cw = run(H, LLRs, max_num_iterations)
    if(decoded_cw[1]):
        if(not np.all((cw==decoded_cw[0].T[0]))):
            undet_err_count = 1
            undet_err_cws.append((cw,decoded_cw[0]))
    else:
        err_count = 1
    cws_all.append((cw,decoded_cw[0]))
    return(bit_err_prob,err_count,undet_err_count,decoded_cw[2],undet_err_cws,cws_all,LLRs)


def LDPC_test_run(H, G_t, bit_err_prob, max_num_iterations, num_trials):
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

##def H_Galg_construct(n, wc, wr):
##    k = n*wc/wr
##    H_sub = np.zeros((n/wr,n))
##    for i in range(1, n/wr+1):
##        for j in range((i-1)*wr+1,i*wr+1):
##            H_sub[i,j] = H_sub[i,j] + 1
##    H_pre = H_sub
##    for  i in range(2,wc+1):
##
##    H = np.zeros((k,n))
##    for i in range(1,wc+1):
##        H

def make_H(n,k, ratio_of_ones):
    I_nk = np.identity(n-k)
    #minus_P_trans = bernoulli.rvs(ratio_of_ones, size=(n-k,k))
    return np.concatenate((minus_P_trans, I_nk),1)


## Obsolete
def make_H_gallager(N,col_w,row_w):
    """
    make_H_gallager(N,col_w,row_w):

    N:     Length of code.
    col_w: Column weight of code (i.e., how many 1's in a column).
    row_w: Row weight of code (i.e., how many 1's in a row).
    
    Create a regular LDPC code matrix using the construction in
    Gallager's book and return the result as a link array.  The
    return value, r, is a list of N lists where r[i][j] is the
    ith one in the jth row of H.

    The algorithm works by starting with a sub-matrix of n/row_w rows
    containing col_w 1's in the first positions in the first row,
    col_w 1's in the next set of positions in the second row, etc.
    
    H = [ 1 1 1 1 0 0 0 0 0 0 0 0 ...
          0 0 0 0 1 1 1 1 0 0 0 0 ...
          0 0 0 0 0 0 0 0 1 1 1 1 ...
          ...........................

    Next we create col_w-1 additional submatricies which are permuations of
    this and append them to the bottom of H.
    """
    num_rows = (N * col_w) / row_w
    num_sub_rows = num_rows / col_w
    assert row_w*num_rows == N*col_w, 'N*col_w not divisible by row_w'
    assert (N/row_w)*row_w == N, 'N not divisible by row_w'

    H_sub = [0]*num_sub_rows
    for i in range(num_sub_rows):
        H_sub[i] = map(lambda x,y: x + y,[i*row_w]*row_w,range(row_w))
    H = copy.deepcopy(H_sub)

    for i in range(col_w-1):

        H_new_sub = [0]*num_sub_rows
        for m in range(num_sub_rows):
            H_new_sub[m] = [0]*row_w

        
        rand_perm = randperm(N)
        for j in range(num_sub_rows):
            for k in range(row_w):
                H_new_sub[j][k] = rand_perm[H_sub[j][k]]
        l = list(H_new_sub[j])
        l.sort()
        H_new_sub[j] = l
        H = H + copy.deepcopy(H_new_sub)


    result = np.array([[0 for j in range(N)] for i in range(num_rows)])
    for i in range(num_rows):
        for j in H[i]:
            result[i][j] = 1

    return result


def randperm(N):
    "randperm(N): Return a random permuation of [0,1,...,N-1]."
    result = [None]*N
    for i in range(N):
        index = random.randrange(N)
        while (result[index] != None):
            index = random.randrange(N)
        result[index] = i
    return result



                

def neighbors_vr(H, vr):
    neighbors = []
    for i in range(len(H.T[vr])):
        if(H.T[vr][i]==1):
            neighbors.append(i)
    return neighbors

def neighbors_cn(H, cn):
    neighbors = []
    for i in range(len(H[cn])):
        if(H[cn][i]==1):
            neighbors.append(i)
    return neighbors

def encode(data, H):
    I_k = np.identity(len(H[0])-len(H))
    P = np.split(H, [len(H[0])-len(H)], 1)[0].T
    G = np.concatenate((I_k, P),1)
    print G
    result_array = np.dot(np.array([int(i) for i in data]),G)%2
    result_string = ''
    for i in result_array:
        result_string += str(int(i))
    return (result_array, result_string)

def LLR_to_cw(LLRs):
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

def update_vr_beliefs(H, cn_to_vr_msgs, LLR_init):
    LLRs = []
    for i in range(len(H[0])):
        belief = LLR_init[i]
        #print(type(belief))
        for j in neighbors_vr(H, i):
            #print(j)
            #print(cn_to_vr_msgs[j][i])
            #print(type(cn_to_vr_msgs[j][i]))
            belief += cn_to_vr_msgs[j][i]
        LLRs.append(belief)
    return LLRs

def update_vr_to_cn_msgs(H, cn_to_vr_msgs, vr_to_cn_msgs, LLR_init):
    for i in range(len(H[0])):
        msg = LLR_init[i]
        for j in neighbors_vr(H, i):
            msg += cn_to_vr_msgs[j][i]
        for j in neighbors_vr(H, i):
            vr_to_cn_msgs[i][j] = msg - cn_to_vr_msgs[j][i]

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
    debug_b=[]
    debug_ch=[]
    debug_bf=[]
    #print('ck1')
    vr_to_cn_msgs = vr_to_cn_msgs_init(LLR_init, H)
    #msg_info = vr_to_cn_msgs
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
        #print LLR_to_cw(LLRs)[0].T
        debug_b.append(deepcopy(vr_to_cn_msgs))
        debug_ch.append(deepcopy(cn_to_vr_msgs))
        debug_bf.append(deepcopy(LLRs))
        if(np.array_equal(np.dot(H,LLR_to_cw(LLRs)[0])%2, np.array([[0] for j in range(len(H))]))):
            #print('ck7')
            return (LLR_to_cw(LLRs)[0]%2, i+1, clock()-t1,debug_b,debug_ch,debug_bf)
        update_vr_to_cn_msgs(H, cn_to_vr_msgs, vr_to_cn_msgs, LLR_init)
        #print('ck8')
    return (LLR_to_cw(LLRs)[0]%2, False, clock()-t1, debug_b,debug_ch,debug_bf)

def check_dec_eq_enc(enc, dec):
    for i in range(len(enc)):
        if(enc[i] != dec[i][0]):
            return (False,i)
    return (True,0)

def test_decoder(bit_err_prob, number_times, n, k, ratio_of_ones):
    if(bit_err_prob==0):
        bit_err_prob = 1e-10
    elif(bit_err_prob >= 1):
        bit_err_prob = 1-(1e-10)
    parity_check_matrix = make_H(n, k, ratio_of_ones)
    CW = ''
    for i in range(k):
        CW += str(random.randint(0,1))
    print('CW = ' + CW)
    codeword = encode(CW, parity_check_matrix)
    
    sector_failure_times = 0

    for i in range(number_times):
        LLRs = []
        noise_word = bernoulli.rvs(bit_err_prob, size=n)
        modulated_cw = (codeword[0] + noise_word)%2
        
        for j in modulated_cw:
            if(j==0):
                LLRs.append(math.log((1-bit_err_prob)/bit_err_prob))
            if(j==1):
                LLRs.append(math.log(bit_err_prob/(1-bit_err_prob)))
        print('H = ' + str(parity_check_matrix))
        print('LLRs = '+str(LLRs))
        print('noise_word = '+str(noise_word))
        decoded_cw = run(parity_check_matrix, np.array(LLRs), 10)
##        print('decoded_cw = ' + str(decoded_cw[0]))
##        print(decoded_cw[1])
##        print(check_dec_eq_enc(codeword[0],decoded_cw[0])[0])
        if((not decoded_cw[1]) or (not check_dec_eq_enc(codeword[0],decoded_cw[0])[0])):
            sector_failure_times += 1.0

    return (sector_failure_times/number_times, parity_check_matrix)

def LDPC_decoder_curve(SNR_range, number_times, n, k, ratio_of_ones):
    SNR_axis = np.linspace(SNR_range[0],SNR_range[1]-SNR_range[2],(SNR_range[1]-SNR_range[0])/SNR_range[2])
    y_axis = [test_decoder(SNR_axis[int(((SNR_range[1]-SNR_range[0])/SNR_range[2])-1-i)], number_times, n, k, ratio_of_ones)[0] for i in range(int((SNR_range[1]-SNR_range[0])/SNR_range[2]))]
    pylab.plot(SNR_axis, y_axis)
    pylab.grid(True)
    pylab.show()
    return [SNR_axis, y_axis]

# John Fan's book example
#decoded_cw = run(H, np.array([-2,-2,-1,2,-2,0]), 10)[0]

##H = make_H(100, 90, 0.2)
#LLRs = []

##CW = ''
##for i in range(90):
##    CW += str(random.randint(0,1))
##print('CW = ' + CW)
##
##codeword = encode(CW, H)
##decoded_cw = run(H, np.array(LLRs), 10)[0]
##

H = parse_MacKay_dense('204.33.484 (N=204,K=102,M=102,R=0.5).txt',102,204)
G_t = np.loadtxt('204.33.484 (N=204,K=102,M=102,R=0.5)G_t.txt')


##if __name__ == '__main__':
##    p = Pool(4)
##    print(p.map(randperm, [10,17, 7, 9]))
    
