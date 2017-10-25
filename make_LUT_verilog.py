import math
from deModel import DeFixedInt as fix


int_width = 7
frac_width = 4

LUT_string = ''

fil = open('tanh_LUT.txt','w')


for i in range(-2**10,2**10):
    LUT_string += 'when ' + str(i) + ' =>\noutput <= "' + fix(int_width,frac_width,math.tanh(i/(2.0**frac_width))).bit() + '"\n'


fil.write(LUT_string)
