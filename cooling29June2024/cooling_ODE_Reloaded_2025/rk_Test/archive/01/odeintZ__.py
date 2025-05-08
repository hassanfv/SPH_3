
import numpy as np


kmax = 100 # We want to save 100 intermediate data.
kount = 0
nvar = len(ystart)

xp = np.zeros(kmax)         # Will contain data for later tests and checks!
yp = np.zeros((nvar, kmax)) # Will contain data for later tests and checks!


eps = None   #!!!!!!!!!!!!!!!!!!! PLACE HOLDER !!!!!!!!!!!!!!
hmin = None  #!!!!!!!!!!!!!!!!!!! PLACE HOLDER !!!!!!!!!!!!!!


def odeint(ystart, x1, x2, h1):

  MAXSTP=10000
  TINY=1.0e-30

  yscal = np.zeros(yscal)
  y = np.zeros(yscal)
  dydx = np.zeros(yscal)

  x = x1

  h = np.sign(h1) * (x2 - x1)

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
        xp[kout] = x
        xsav = x
        kount += 1

    if (x+h-x2)*(x+h-x1) > 0.0:
      h = x2 - x

    x, y, hdid, hnext = rkqs(y,dydx,x,h,eps,yscal,hdid,hnext,derivs)

    if hdid == h:
      nok += 1
    else
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
    
    if abs(hnext) <= hmin:
      raise ValueError("Step size too small in odeint")

  raise ValueError("Too many steps in routine odeint")

  





