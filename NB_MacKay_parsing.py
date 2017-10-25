fil = open('qexp.txt')
strin = fil.read()
strin_l = strin.split('\n')
strin_new = []

for i in strin_l:
    i += '0 0'
    strin_new.append(i)

strin_new_l = '\n'.join(strin_new)
new_fil = open('qw.txt','w')
new_fil.write(strin_new_l)
