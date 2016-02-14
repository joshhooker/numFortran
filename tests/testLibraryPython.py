import h5py
import numpy as np
import scipy.special as sp

def drange(start, stop, step):
  r = start
  while r < stop:
    yield r
    r += step

# Bessel Function of the 1st Kind
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

# Bessel Function of the 2nd Kind
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

# Modified Bessel Function of the 1st Kind
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

# Modified Bessel Function of the 2nd Kind
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

# Spherical Bessel Function of the 1st Kind
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

# Spherical Bessel Function of the 2nd Kind
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

# Airy A Function
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

# Airy B Function
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

# Error Function
err_x=[]
err=[]
for x in drange(-5.0,5.0,0.2):
  err_x.append(x)
  err.append(sp.erf(x))
err_length = len(err)
f = open('test_err.dat','w')
f.write('{}\n'.format(err_length))
for i in range (0,err_length):
  f.write('{} {}\n'.format(err_x[i],err[i]))
f.close()

# Complementary Error Function
errc_x=[]
errc=[]
for x in drange(-5.0,5.0,0.2):
  errc_x.append(x)
  errc.append(sp.erfc(x))
errc_length = len(err)
f = open('test_errc.dat','w')
f.write('{}\n'.format(errc_length))
for i in range (0,errc_length):
  f.write('{} {}\n'.format(errc_x[i],errc[i]))
f.close()

# Beta Function
beta_a=[]
beta_b=[]
beta=[]
for a in drange(0.1,5.0,0.1):
  for b in drange(0.1,5.0,0.1):
    beta_a.append(a)
    beta_b.append(b)
    beta.append(sp.beta(a,b))
beta_length = len(beta)
f = open('test_beta.dat','w')
f.write('{}\n'.format(beta_length))
for i in range (0,beta_length):
  f.write('{} {} {}\n'.format(beta_a[i],beta_b[i],beta[i]))
f.close()

