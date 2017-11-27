import numpy as np

f = open("noise.txt")

x = f.readline()
x = x.split(" ")
x = x[:-1]
x = np.array(x)
y = f.readline()
y = y.split(" ")
y = y[:-1]
y = np.array(y)

indices = []
for i in range(13):
    indices.append(np.where(x==str(i))[0])
    
