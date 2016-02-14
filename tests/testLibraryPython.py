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

in_n=[]
in_x=[]
in_arr=[]
for n in drange(-3.0,5.0,0.2):
  for x in drange(1.0,5.0,0.2):
    in_n.append(n)
    in_x.append(x)
    in_arr.append(sp.iv(n,x))
in_length = len(in_arr)
f = open('test_in.dat','w')
f.write('{}\n'.format(in_length))
for i in range(0,in_length):
  f.write('{} {} {}\n'.format(in_n[i],in_x[i],in_arr[i]))
f.close()

kn_n=[]
kn_x=[]
kn=[]
for n in drange(-3.0,5.0,0.2):
  for x in drange(1.0,5.0,0.2):
    kn_n.append(n)
    kn_x.append(x)
    kn.append(sp.kv(n,x))
kn_length = len(kn)
f = open('test_kn.dat','w')
f.write('{}\n'.format(kn_length))
for i in range(0,kn_length):
  f.write('{} {} {}\n'.format(kn_n[i],kn_x[i],kn[i]))
f.close()

sphjn_n=[]
sphjn_x=[]
sphjn=[]
for n in range(0,10,1):
  for x in drange(0.0,10.0,0.2):
    sphjn_n.append(n)
    sphjn_x.append(x)
    sphjn.append(sp.sph_jn(n,x)[0][n])
sphjn_length = len(sphjn)
f = open('test_sphjn.dat','w')
f.write('{}\n'.format(sphjn_length))
for i in range(0,sphjn_length):
  f.write('{} {} {}\n'.format(sphjn_n[i],sphjn_x[i],sphjn[i]))
f.close()

sphyn_n=[]
sphyn_x=[]
sphyn=[]
for n in range(0,10,1):
  for x in drange(0.2*n+0.2,10.0,0.2):
    sphyn_n.append(n)
    sphyn_x.append(x)
    sphyn.append(sp.sph_yn(n,x)[0][n])
sphyn_length = len(sphyn)
f = open('test_sphyn.dat','w')
f.write('{}\n'.format(sphyn_length))
for i in range(0,sphyn_length):
  f.write('{} {} {}\n'.format(sphyn_n[i],sphyn_x[i],sphyn[i]))
f.close()

airyA_x=[]
airyA=[]
for x in drange(-5.0,5.0,0.2):
  airyA_x.append(x)
  airyA.append(sp.airy(x)[0])
airyA_length = len(airyA)
f = open('test_airyA.dat','w')
f.write('{}\n'.format(airyA_length))
for i in range (0,airyA_length):
  f.write('{} {}\n'.format(airyA_x[i],airyA[i]))
f.close()

airyB_x=[]
airyB=[]
for x in drange(-5.0,5.0,0.2):
  airyB_x.append(x)
  airyB.append(sp.airy(x)[2])
airyB_length = len(airyB)
f = open('test_airyB.dat','w')
f.write('{}\n'.format(airyB_length))
for i in range (0,airyB_length):
  f.write('{} {}\n'.format(airyB_x[i],airyB[i]))
f.close()
