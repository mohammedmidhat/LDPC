import numpy as np

grey_code = [15, 14, 10, 8, 9, 1, 0, 2, 6, 4, 12, 13, 5, 7, 3, 11]

num_strs = 64*4
num_bytes_per_page = 18336
num_pages_per_str = 4

f = open('data.bin','wb')

rand_data = np.random.randint(2, 15, size=8*num_bytes_per_page*num_strs)

lp_data = [0 for i in range(num_bytes_per_page)]
mp_data = [0 for i in range(num_bytes_per_page)]
up_data = [0 for i in range(num_bytes_per_page)]
tp_data = [0 for i in range(num_bytes_per_page)]

for k in range(num_strs):
    for i in range(num_bytes_per_page):
        lp_byte = 0
        mp_byte = 0
        up_byte = 0
        tp_byte = 0
        
        for j in range(8):
            symbol = rand_data[k*num_bytes_per_page*8+i*8+j]
            symbol_code = grey_code[symbol]
            
            lp_bit = symbol_code & 1
            mp_bit = (symbol_code >> 1) & 1
            up_bit = (symbol_code >> 2) & 1
            tp_bit = (symbol_code >> 3) & 1

            lp_byte += lp_bit << (8-1-j)
            mp_byte += mp_bit << (8-1-j)
            up_byte += up_bit << (8-1-j)
            tp_byte += tp_bit << (8-1-j)

        lp_data[i] = lp_byte
        mp_data[i] = mp_byte
        up_data[i] = up_byte
        tp_data[i] = tp_byte

    f.write(bytearray(lp_data))
    f.write(bytearray(mp_data))
    f.write(bytearray(up_data))
    f.write(bytearray(tp_data))
    
f.close()

##for i in range(num_strs):
##    str_data = file_data[i*num_pages_per_str*num_bytes_per_page:(i+1)*num_pages_per_str*num_bytes_per_page]
##    for j in range(num_bytes_per_page):
##        lp_byte = str_data[j]
##        mp_byte = str_data[j+num_bytes_per_page]
##        up_byte = str_data[j+2*num_bytes_per_page]
##        tp_byte = str_data[j+3*num_bytes_per_page]
##
##        for k in range(8):
##            (8-1-i)
