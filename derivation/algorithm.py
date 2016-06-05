# This repo represents a closer-to-viable rearrangement of the worked example
# in derivation.py

import numpy as np
from math import factorial

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

  brn = bn / n
  comp1 = np.sum(brn * np.power(rtest1 / r, np.arange(n)))
  comp2 = np.sum(brn * np.power(rtest2 / r, np.arange(n)))
  comp3 = np.sum(brn * np.power(rtest3 / r, np.arange(n)))
  print ftest1
  print comp1

  diff1 = comp1 - ftest1
  diff2 = comp2 - ftest2
  diff3 = comp3 - ftest3

  max_error = np.max(np.abs([diff1, diff2, diff3])) / np.max(np.abs([ftest1, ftest2, ftest3]))

  return max_error > 1e-3

def circle(z, r, n):
  return z + r * np.exp(np.linspace(0.0, np.pi * 2.0j, num=n, endpoint=False))

cnt = 0
def f(z):
  global n
  global cnt
  cnt = cnt + n
  #return np.exp(1.0j * z)
  #return np.tan(z)
  #return 1.0j + z + 1.0j * z**2
  return 1.0 / (1.0 - z)

N = 1
direction_changes = 0
rs = []
bs = []
z0 = 0.0
r = 0.6580924658
fac = 2
pdirec = None
degenerate = False

if N <= 6:
  range_num = 1
  n = 8
elif N <= 12:
  range_num = 2
  n = 16
elif N <= 25:
  range_num = 3
  n = 32
elif N <= 51:
  range_num = 4
  n = 64
else:
  range_num = 5
  n = 128

crat = np.exp(np.log(1e-4) / (n - 1))

degenerate = False

for iters in xrange(100):
  print 'r = %g' % (r)
  iters += 1

  bn = np.fft.fft(f(circle(z0, r, n)))

  if direction_changes > 1 or degenerate:
    bs.append(bn / n * np.power(r, -np.arange(n)))
    rs.append(r)
    if len(rs) >= 3:
      break

  if not degenerate:
    bnc = bn / crat ** np.arange(n)
    m1 = np.max(np.abs(bnc[:n // 2]))
    m2 = np.max(np.abs(bnc[n // 2:]))

    degenerate =  m1 / m2 < 1e-8 or m2 / m1 < 1e-8

  if degenerate:
    needs_smaller = iters % 2 == 0
  else:
    needs_smaller = m1 < m2 or test2(z0, r, f, n, bn)

  if pdirec is not None and needs_smaller != pdirec:
    direction_changes += 1

  if direction_changes > 0:
    fac = np.sqrt(fac)

  if needs_smaller:
    r /= fac
  else:
    r *= fac

  pdirec = needs_smaller

print np.real(bs[0])


extrap1 = bs[1] - (bs[1] - bs[0]) / (1.0 - (rs[0] / rs[1])**n)
extrap2 = bs[2] - (bs[2] - bs[1]) / (1.0 - (rs[1] / rs[2])**n)
extrap3 = extrap2 - (extrap2 - extrap1) / (1.0 - (rs[0] / rs[2])**n)

extrap3 *= [factorial(i) for i in range(len(extrap3))]

print 'Function evaluations: ', cnt

print 'answer:'
for i in range(len(extrap3)):
  print '%3i: %24.18f + %24.18fj' % (i, np.real(extrap3[i]), np.imag(extrap3[i]))
