import numpy as np
import scipy.special as sp

def drange(start, stop, step):
  r = start
  while r < stop:
    yield r
    r += step

f = open('test_jn.dat','w')
for n in drange(-3.0,5.0,0.2):
  for x in drange(1.0,5.0,0.2):
    f.write('{} {} {}\n'.format(n,x,sp.jv(n,x)))
f.close()

f = open('test_yn.dat','w')
for n in drange(-3.0,5.0,0.2):
  for x in drange(1.0,5.0,0.2):
    f.write('{} {} {}\n'.format(n,x,sp.yv(n,x)))
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
  for x in drange(0.2*n,10.0,0.2):
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