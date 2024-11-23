import numpy as np
import matplotlib.pyplot as plt
import time

#----- rk4
def rk4(y, dydx, x, h):

  n = len(y)
  
  yout = np.zeros(n)
  dyt = np.zeros(n)
  yt = np.zeros(n)
  dym = np.zeros(n)
  
  hh = 0.5 * h
  h6 = h / 6.0
  
  xh = x + hh # i.e. x + h/2
  
  dydx = derivs(x, y) # this is k1
  
  for i in range(n):
    yt[i] = y[i] + hh * dydx[i] # y0 + h/2 * k1
  dyt = derivs(xh, yt) # this is k2.
  
  for i in range(n):
    yt[i] = y[i] + hh * dyt[i]
  dym = derivs(xh, yt) # this is k3.
  
  for i in range(n):
    yt[i] = y[i] + h * dym[i]
    dym[i] += dyt[i] # clever idea !
  dyt = derivs(x + h, yt)
  
  for i in range(n):
    yout[i] = y[i] + h6 * (dydx[i] + dyt[i] + 2.0 * dym[i])
  
  return yout
  


def derivs(x, y):
  dydx = np.zeros_like(y)
  dydx[0] = -0.013*y[0]-1000.0*y[0]*y[2];
  dydx[1] = -2500.0*y[1]*y[2];
  dydx[2] = -0.013*y[0]-1000.0*y[0]*y[2]-2500.0*y[1]*y[2];
  return dydx

h = 0.001

x = 0.0
y = np.array([1e-4, 1e-3, 1e-3])

dydx = derivs(x, y)

TA = time.time()
a = rk4(y, dydx, x, h)
print('Elapsed time = ', time.time() - TA)

print(a)




