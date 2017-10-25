import numpy as np

def matrix_multip_init(n,k):
    result = 'module matrix_multip #(parameter WIDTH = 20, parameter N = ' + str(n) + ', parameter K = ' + str(k) + ')(\n'
    result += 'input [N*(N-K)-1:0] H_t,\n'
    result += 'input [N-1:0] cw,\n'
    result += 'output dot_product\n'
    result += ');\n'
    result += '\n'

    for i in range(n):
        result += 'wire [N-K-1:0] element_' + str(i) + ';\n'
        result += 'assign element_' + str(i) + ' = H_t[' + str(i+1) + '*(N-K)-1:' + str(i) + '*(N-K)]*cw[' + str(i) + '];\n'

    result += 'wire [N-K-1:0] dot_product_vector;\n'

    result += 'assign dot_product_vector = '
    for i in range(n):
        result += 'element_' + str(i) + '^'
    result = result[:-1] + ';\n'
    result += 'assign dot_product = (dot_product_vector)? 1:0;\n'

    result += 'endmodule\n'

    return result


def min_init(cn_deg):
    result = 'module min #(parameter WIDTH = 11)(\n'

    for i in range(1,cn_deg):
        result += 'input [WIDTH-1:0] msg_' + str(i) + ',\n'

    result += 'output [WIDTH-1:0] msg\n'
    result += ');\n'
    result += '\n'

    for i in range(1,cn_deg):
        result += 'wire [WIDTH-1:0] abs_msg_' + str(i) + ';\n'
        result += 'assign abs_msg_' + str(i) + ' = (msg_' + str(i) + '[WIDTH-1])? -msg_' + str(i) + ':msg_' + str(i) + ';\n'
    result += 'wire [WIDTH-1:0] abs_msg;\n'


def LLRs_to_cw_init(n):
    result = 'module LLRs_to_cw #(parameter WIDTH = 11, parameter N = ' + str(n) + ')(\n'
    result += 'input [N*WIDTH-1:0] LLRs,\n'
    result += 'output [N-1:0] result\n'
    result += ');\n'
    result += '\n'

    for i in range(n):
        result += 'assign result[' + str(i) + '] = (LLRs[' + str(i+1) + '*WIDTH-1])? 1:0;\n'

    result += 'endmodule\n'

    return result

def vn_init(vn_deg):
    result = 'module vn #(parameter WIDTH = 11)(\n'

    for i in range(vn_deg):
        result += 'input [WIDTH-1:0] msg_in_' + str(i+1) + ',\n'
    for i in range(vn_deg):
        result += 'output [WIDTH-1:0] msg_out_' + str(i+1) + ',\n'

    result += 'input [WIDTH-1:0] belief\n'
    result += ');\n'
    result += '\n'

    for i in range(vn_deg):
        result += 'assign msg_out_' + str(i+1) + ' = belief - msg_in_' + str(i+1) + ';\n'

    result += 'endmodule\n'

    return result


def cn_init(cn_deg):
    result = 'module cn #(parameter WIDTH = 11)(\n'

    for i in range(cn_deg):
        result += 'input [WIDTH-1:0] msg_in_' + str(i+1) + ',\n'
    for i in range(cn_deg):
        result += 'output [WIDTH-1:0] msg_out_' + str(i+1) + ',\n'

    result = result[:-3] + '\n);\n'
    result += '\n'

    for i in range(cn_deg):
        result += 'min min_to_v' + str(i+1) + '_inst (\n'
        for j in range(1,cn_deg):
            result += '.msg_' + str(j) + ' (msg_in_' + str(j if(j<i) else j+1) + '),\n'
        result += '.msg (msg_out_' + str(i+1) + ')\n'
        result += ');\n'
        result += '\n'

    result += 'endmodule\n'

    return result


def posterior_belief_init(vn_deg):
    result = 'module posterior_belief #(parameter WIDTH = 11)(\n'

    for i in range(vn_deg):
        result += 'input [WIDTH-1:0] msg_in_' + str(i+1) + ',\n'

    result += 'input [WIDTH-1:0] LLR_init,\n'
    result += 'output [WIDTH-1:0] result\n'
    result += ');\n'
    result += '\n'

    result += 'assign result = LLR_init'
    for i in range(vn_deg):
        result += ' + msg_in_' + str(i+1)
    result += ';\n'

    result += 'endmodule\n'

    return result

    
def vns_init(H, n, c_w):
    result = '///////////////////////////////////////////////////////\n'
    result += '// Variable Nodes\n'
    result += '///////////////////////////////////////////////////////\n'
    for i in range(n):
        result += 'vn vn_' + str(i) + '_inst (\n    '
        for j in range(c_w):
            result += '.msg_in_' + str(j+1) + '(cn_vr_msgs[' + str(H[i,j]) + '][' + str(i) + ']),\n    '
        result += '.belief(' +  + '),\n    '
        for j in range(c_w):
            if(j==c_w-1):
                result += '.msg_out_' + str(j+1) + '(vn_' + str(i) + '_out_' + str(j+1) + ')\n    '
            else:
                result += '.msg_out_' + str(j+1) + '(vn_' + str(i) + '_out_' + str(j+1) + '),\n    '
        result += ');\n\n'
    return result


def posterior_beliefs_init(H, n, c_w):

    result = '///////////////////////////////////////////////////////\n'
    result += '// Posterior Belief\n'
    result += '///////////////////////////////////////////////////////\n'

    for i in range(n):
        result += 'posterior_belief posterior_belief_' + str(i) + '_inst (\n    '
        for j in range(c_w):
            result += '.msg_in_' + str(i+1) + '(cn_vr_msgs[' + str(H[i,j]) + '][' + str(i) + ']),\n    '
        result += '.LLR_init(LLRs[' + str(i) + ']),\n    '

    
def LLRs_to_cw_init(n):
    result = 'module LLRs_to_cw #(parameter WIDTH = 11, parameter N = ' + str(n) + ')(\n'
    result += 'input [N*WIDTH-1:0] LLRs,\n'
    result += 'output [N-1:0] result\n'
    result += ');\n'
    result += '\n'

    for i in range(n):
        result += 'assign result[' + str(i) + '] = (LLRs[' + str(i+1) + '*WIDTH-1])? 1:0;\n'

    result += 'endmodule\n'

    return result


fi = open('exp_fil.txt','w')
fi.write(matrix_multip_init(204,102))
fi.close()
