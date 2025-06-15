
import numpy as np
import matplotlib.pyplot as plt
import time


#----- sign
def sign(a, b):
  if b >= 0:
    return a if a >= 0 else -a
  else:
    return -a if a >= 0 else a


#===== rkck
def rkck(y, dydx, x, h):
  a2, a3, a4, a5, a6 = 0.2, 0.3, 0.6, 1.0, 0.875
  b21 = 0.2
  b31, b32 = 3.0 / 40.0, 9.0 / 40.0
  b41, b42, b43 = 0.3, -0.9, 1.2
  b51, b52, b53, b54 = -11.0 / 54.0, 2.5, -70.0 / 27.0, 35.0 / 27.0
  b61, b62, b63, b64, b65 = 1631.0 / 55296.0, 175.0 / 512.0, 575.0 / 13824.0, 44275.0 / 110592.0, 253.0 / 4096.0
  c1, c3, c4, c6 = 37.0 / 378.0, 250.0 / 621.0, 125.0 / 594.0, 512.0 / 1771.0
  dc1 = c1 - 2825.0 / 27648.0
  dc3 = c3 - 18575.0 / 48384.0
  dc4 = c4 - 13525.0 / 55296.0
  dc5 = -277.0 / 14336.0
  dc6 = c6 - 0.25

  n = len(y)
  ak2 = np.zeros(n)
  ak3 = np.zeros(n)
  ak4 = np.zeros(n)
  ak5 = np.zeros(n)
  ak6 = np.zeros(n)
  ytemp = np.zeros(n)
  yout = np.zeros(n)
  yerr = np.zeros(n)

  for i in range(n):
      ytemp[i] = y[i] + b21 * h * dydx[i]
  ak2 = derivs(x + a2 * h, ytemp)

  for i in range(n):
      ytemp[i] = y[i] + h * (b31 * dydx[i] + b32 * ak2[i])
  ak3 = derivs(x + a3 * h, ytemp)

  for i in range(n):
      ytemp[i] = y[i] + h * (b41 * dydx[i] + b42 * ak2[i] + b43 * ak3[i])
  ak4 = derivs(x + a4 * h, ytemp)

  for i in range(n):
      ytemp[i] = y[i] + h * (b51 * dydx[i] + b52 * ak2[i] + b53 * ak3[i] + b54 * ak4[i])
  ak5 = derivs(x + a5 * h, ytemp)

  for i in range(n):
      ytemp[i] = y[i] + h * (b61 * dydx[i] + b62 * ak2[i] + b63 * ak3[i] + b64 * ak4[i] + b65 * ak5[i])
  ak6 = derivs(x + a6 * h, ytemp)

  for i in range(n):
      yout[i] = y[i] + h * (c1 * dydx[i] + c3 * ak3[i] + c4 * ak4[i] + c6 * ak6[i])

  for i in range(n):
      yerr[i] = h * (dc1 * dydx[i] + dc3 * ak3[i] + dc4 * ak4[i] + dc5 * ak5[i] + dc6 * ak6[i])

  return yout, yerr


#===== rkqs
def rkqs(y, dydx, x, htry, yscal):

  SAFETY = 0.9
  PGROW = -0.2
  PSHRNK = -0.25
  ERRCON = 1.89e-4

  n = len(y)
  
  h = htry
  
  while True:
    ytemp, yerr = rkck(y, dydx, x, h)
    
    errmax = 0.0
    for i in range(n):
      errmax = max(errmax, abs(yerr[i] / yscal[i]))

    errmax /= eps
    if errmax <= 1.0:  # Accept h if error is within tolerance.
      break
    
    htemp = SAFETY * h * errmax ** PSHRNK # If h is not within tolerance, we try to decrease h.
    h = max(htemp, 0.1 * h) if h >= 0 else min(htemp, 0.1 * h) # making sure that htemp is never below 0.1 times the previous h.
    
    xnew = x + h # I think xnew is only used to check whether "xnew = x" (meaning h is 0) so that we raise the "stepsize underflow" error !
    
    if xnew == x:  # Stepsize underflow ---> This is used to alert us if h = 0.0 !
      raise ValueError("Stepsize underflow in rkqs")

  # When we exit the WHILE loop, we have the "y" (ytemp) which is evolved for one step h with the best h (i.e. with a h that gives errmax <= 1.0)!
  if errmax > ERRCON:
    hnext = SAFETY * h * errmax ** PGROW # Here we decrease the stepsize.
  else:
    hnext = 5.0 * h # Here we increase the stepsize.

  x += h
  hdid = h

  for i in range(n):
    y[i] = ytemp[i]

  return x, y, hdid, hnext


#===== odeint
def odeint(ystart, x1, x2, h1):

  MAXSTP=10000
  TINY=1.0e-30

  kount = 0

  nvar = len(ystart)

  yscal = np.zeros(nvar)
  y = np.zeros(nvar)
  dydx = np.zeros(nvar)

  x = x1

  #h = np.sign(h1) * (x2 - x1) # This was wrong. I made a worng conversion from C++ to python !
  h = sign(h1, x2 - x1)

  nok = nbad = 0

  for i in range(nvar):
    y[i] = ystart[i]

  if kmax > 0:
    xsav = x - dxsav*2.0

  for nstp in range(MAXSTP): # main for loop !
    dydx = derivs(x, y)
    for i in range(nvar):
      yscal[i] = abs(y[i]) + abs(dydx[i]*h) + TINY;

    if kmax > 0 and kount < kmax - 1 and abs(x-xsav) > abs(dxsav):
      for i in range(nvar):
        yp[i, kount] = y[i]
        xp[kount] = x
        xsav = x
      kount += 1

    if (x+h-x2)*(x+h-x1) > 0.0: # To adjust h if x+h exceed the upper limit value, x2 !
      h = x2 - x

    x, y, hdid, hnext = rkqs(y, dydx, x, h, yscal)

    if hdid == h:
      nok += 1
    else:
      nbad += 1

    if (x-x2)*(x2-x1) >= 0.0:
      for i in range(nvar):
        ystart[i] = y[i]
      if kmax != 0:
        for i in range(nvar):
          yp[i, kount] = y[i]
        xp[kount] = x
        kount += 1
      return ystart, xp, yp, nok, nbad
    
    print('x, hnext, hmin = ', x, hnext, hmin)
    
    if abs(hnext) <= hmin:
      raise ValueError("Step size too small in odeint")

  raise ValueError("Too many steps in routine odeint")


# ===== derivs (Lorenz butter fly feature)
def derivs(x, y):
    sigma = 10.0
    rho = 28.0
    beta = 8.0 / 3.0

    dydx = np.zeros_like(y)
    dydx[0] = sigma * (y[1] - y[0])              # dx/dt = sigma*(y - x)
    dydx[1] = y[0] * (rho - y[2]) - y[1]         # dy/dt = x*(rho - z) - y
    dydx[2] = y[0] * y[1] - beta * y[2]          # dz/dt = x*y - beta*z
    return dydx



kmax = 2000 # We want to save 100 intermediate data. !!!!!!! Adjust if need or less sampling points !!!!!!!!
dxsav = 0.001 #!!!!!!!!!!!!!!!!!!!!!!!!!! To be adjusted for each problem !!!!!!!!!!!!!!!!!!!!!!!!!!


# Note that here x1 and x2 represent time !
x1 = 0.0 # Initial time !!!!!!!!!!!!!!!!!!!!!!!
x2 = 20.0 # Final time !!!!!!!!!!!!!!!!!!!!!!!!!
h1 = 0.1 # Step size !!!!!!!!!!!!!!!!!!!!!!!!
ystart = np.array([1.0, 1.0, 1.0]) # Initial condition

nvar = len(ystart)

xp = np.zeros(kmax)         # Will contain data for later tests and checks!
yp = np.zeros((nvar, kmax)) # Will contain data for later tests and checks!

eps = 1e-7   # THIS affects the execution time !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
hmin = 1e-5  #!!!!!!!!!!!!!!!!!!!!!!!!!! To be adjusted for each problem !!!!!!!!!!!!!!!!!!!!!!!!!!

TA = time.time()
ystart, xp, yp, nok, nbad = odeint(ystart, x1, x2, h1)
print('Elapsed time = ', time.time() - TA)

#print(yp.shape)
#print(xp)
print()
print('nok, nbad = ', nok, nbad)

y1 = np.zeros(kmax)
y2 = np.zeros(kmax)
y3 = np.zeros(kmax)

for j in range(kmax):
  y1[j] = yp[0, j]
  y2[j] = yp[1, j]
  y3[j] = yp[2, j]

plt.scatter(y1, y3, s = 5, color = 'r')
plt.plot(y1, y3, linewidth = 0.5, color = 'k')

plt.savefig('fig.png', bbox_inches = 'tight', dpi = 200)

plt.show()





