f_1 = open('verilog_viv.txt')
f_2 = open('verilog.txt')

for i in range(145252):
    if(f_1.readline() != f_2.readline()):
        print(i+1)


f_1.close()
f_2.close()
