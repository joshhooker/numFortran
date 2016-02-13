import h5py
import numpy as np
import scipy.special as sp

def drange(start, stop, step):
  r = start
  while r < stop:
    yield r
    r += step

# f = open('test_jn.dat','w')
# for n in drange(-3.0,5.0,0.2):
#   for x in drange(1.0,5.0,0.2):
#     f.write('{} {} {}\n'.format(n,x,sp.jv(n,x)))
# f.close()
jn_n=[]
jn_x=[]
jn=[]
for n in drange(-3.0,5.0,0.2):
  for x in drange(1.0,5.0,0.2):
    jn_n.append(n)
    jn_x.append(x)
    jn.append(sp.jv(n,x))
jn_length = len(jn)
f = open('test_jn.dat','w')
f.write('{}\n'.format(jn_length))
for i in range(0,jn_length):
  f.write('{} {} {}\n'.format(jn_n[i],jn_x[i],jn[i]))
f.close()

yn_n=[]
yn_x=[]
yn=[]
for n in drange(-3.0,5.0,0.2):
  for x in drange(1.0,5.0,0.2):
    yn_n.append(n)
    yn_x.append(x)
    yn.append(sp.yv(n,x))
yn_length = len(yn)
f = open('test_yn.dat','w')
f.write('{}\n'.format(yn_length))
for i in range(0,yn_length):
  f.write('{} {} {}\n'.format(yn_n[i],yn_x[i],yn[i]))
f.close()

f = open('test_in.dat','w')
for n in drange(-3.0,5.0,0.2):
  for x in drange(1.0,5.0,0.2):
    f.write('{} {} {}\n'.format(n,x,sp.iv(n,x)))
f.close()

f = open('test_kn.dat','w')
for n in drange(-3.0,5.0,0.2):
  for x in drange(1.0,5.0,0.2):
    f.write('{} {} {}\n'.format(n,x,sp.kv(n,x)))
f.close()

f = open('test_sphjn.dat','w')
for n in range(0,10,1):
  for x in drange(0.0,10.0,0.2):
    f.write('{} {} {}\n'.format(n,x,sp.sph_jn(n,x)[0][n]))
f.close()

f = open('test_sphyn.dat','w')
for n in range(0,10,1):
  for x in drange(0.2*n+0.2,10.0,0.2):
    f.write('{} {} {}\n'.format(n,x,sp.sph_yn(n,x)[0][n]))
f.close()

f = open('test_airyA.dat','w')
for x in drange(-5.0,5.0,0.2):
  f.write('{} {}\n'.format(x,sp.airy(x)[0]))
  #print x, sp.airy(x)[0], sp.airy(x)[2]
f.close()

f = open('test_airyB.dat','w')
for x in drange(-5.0,5.0,0.2):
  f.write('{} {}\n'.format(x,sp.airy(x)[2]))
  #print x, sp.airy(x)[0], sp.airy(x)[2]
f.close()