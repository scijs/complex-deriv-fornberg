# This script is simply a full reproduction of the worked example in Fornberg's paper,
# B. Fornberg, "Numerical Differentiation of Analytic Functions", (1981).

import numpy as np

def test1(b, n):
  c0 = 1
  cn1 = 1e-4
  crat = np.exp(np.log(cn1 / c0) / n)
  c = c0 * crat ** np.arange(n)
  print 'test1:', np.argmax(np.abs(b) / c)
  return np.argmax(np.abs(b) / c) <= n // 2

def test2(z, r, f, n, bn):
  rtest1 = r * (-0.4 + 0.3j)
  rtest2 = r * (0.7 + 0.2j)
  rtest3 = r * (0.02 - 0.06j)
  ztest1 = z + rtest1
  ztest2 = z + rtest2
  ztest3 = z + rtest3

  ftest1 = f(ztest1)
  ftest2 = f(ztest2)
  ftest3 = f(ztest3)

  brn = bn * np.power(r, -np.arange(n))
  comp1 = np.sum(brn * np.power((ztest1 - z) , np.arange(n)))
  comp2 = np.sum(brn * np.power((ztest2 - z) , np.arange(n)))
  comp3 = np.sum(brn * np.power((ztest3 - z) , np.arange(n)))

  diff1 = comp1 - ftest1
  diff2 = comp2 - ftest2
  diff3 = comp3 - ftest3

  max_error = np.max(np.abs([diff1, diff2, diff3])) / np.max(np.abs([ftest1, ftest2, ftest3]))
  print 'test2:', max_error
  return max_error < 1e-3

def ensmallen(z, r, f, n, bn):
  if test1(bn, n):
    print 'enlargen 1'
    return False

    if test2(z, r, f, n, bn):
      print 'enlargen 2'
      return False

  print 'ensmallen'
  return True

def circle(z, r, n):
  return z + r * np.exp(np.linspace(0.0, np.pi * 2.0j, num=n, endpoint=False))

def f(z):
  return 1.0 / (1.0 - z)

direction_changes = 0
sample_rads = []
sample_vals = []
range_num = 2
n = 32
z0 = 0.0
r0 = 0.6580924658

#while len(samples) < 3:

a0 = f(circle(z0, r0, n))
b0n = np.fft.fft(a0) / n

print b0n
ensmallen(z0, r0, f, n, b0n)

fac1 = 2
r1 = r0 * fac1
a1 = f(circle(z0, r1, n))
b1n = np.fft.fft(a1) / n

print r1
#print b1n
ensmallen(z0, r1, f, n, b1n)

fac2 = np.sqrt(fac1)
r2 = r1 / fac2
a2 = f(circle(z0, r2, n))
b2n = np.fft.fft(a2) / n

print r2
#print b2n
ensmallen(z0, r2, f, n, b2n)

fac3 = np.sqrt(fac2)
r3 = r2 / fac3
a3 = f(circle(z0, r3, n))
b3n = np.fft.fft(a3) / n

print r3
#print b3n
ensmallen(z0, r3, f, n, b3n)

fac4 = np.sqrt(fac3)
r4 = r3 / fac4
a4 = f(circle(z0, r4, n))
b4n = np.fft.fft(a4) / n

print r4
#print b4n
ensmallen(z0, r4, f, n, b4n)

fac5 = np.sqrt(fac4)
r5 = r4 * fac5
a5 = f(circle(z0, r5, n))
b5n = np.fft.fft(a5) / n

print r5
#print b5n
ensmallen(z0, r5, f, n, b5n)

fac6 = np.sqrt(fac5)
r6 = r5 * fac6
a6 = f(circle(z0, r6, n))
b6n = np.fft.fft(a6) / n

print r6
#print b6n
ensmallen(z0, r6, f, n, b6n)

fac7 = np.sqrt(fac6)
r7 = r6 / fac7
a7 = f(circle(z0, r7, n))
b7n = np.fft.fft(a7) / n

print r7
#print b7n


b5r5n = b5n * np.power(r5, -np.arange(n))
b6r6n = b6n * np.power(r6, -np.arange(n))
b7r7n = b7n * np.power(r7, -np.arange(n))

extrap1 = b6r6n - (b6r6n - b5r5n) / (1.0 - (r5 / r6)**n)
extrap2 = b7r7n - (b7r7n - b6r6n) / (1.0 - (r6 / r7)**n)
extrap3 = extrap2 - (extrap2 - extrap1) / (1.0 - (r5 / r7)**n)

for i in range(len(extrap3)):
  print '%.20f' % np.real(extrap3[i])
