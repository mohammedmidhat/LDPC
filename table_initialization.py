import math

Na = 175.0
Ma = 1023.0
Nf = 173.0
Mf = 1013.0

a_init_list = []
f_init_list = []

a2f = []
f2a = []

for p in range(1,13):
    F_init = Nf*math.log((1-p)*1.0/p,2)
    A_init = -Na*math.log(1-2*p,2)
    a_init_list.append(int(A_init))
    f_init_list.append(int(F_init))

for i in range(1,Ma+1):
    a = 1.0*math.exp(-i/Na)
    f = (1+a)/(1-a)
    z = 1.0*math.log(f,2)
    j = Nf*z+0.5
    if (j <= 0):
        j = 1
    if (j > Mf):
        j = Mf
    a2f.append(int(j))

for i in range(1,Mf+1):
    f = 1.0math.exp(i/Nf)
    a = (f-1)/(f+1)
    z = -1.0*math.log(a,2)
    j = Na*z+0.5
    if (j <= 0):
        j = 1
    if (j > Ma):
        j = Ma
    f2a.append(int(j))

