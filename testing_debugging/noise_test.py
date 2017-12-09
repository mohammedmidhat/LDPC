import numpy as np
import math


num_reads = 2
Q = 13
data_len = 8000

f = open("noise.txt")

x = f.readline()
x = x.split(" ")
x = x[:-1]
x = np.array(x)
y = f.readline()
y = y.split(" ")
y = y[:-1]
y = np.array(y)
logfna = []
for i in range(data_len):
    z = f.readline()
    z = z.split(" ")
    z = z[:-1]
    logfna.append(z)

indices = []
for i in range(Q):
    indices.append(np.where(x==str(i))[0])
    
for i in range(13):
    for j in range(num_reads*Q):
        c = 0
        for k in y[indices[i]]:
            if(k == str(j)):
                c += 1
        #print 1.0*c/len(y[indices[i]])
        
