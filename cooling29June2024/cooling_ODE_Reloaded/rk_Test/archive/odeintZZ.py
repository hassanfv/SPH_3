
#----- sign
def sign(a, b):
  if b >= 0:
    return a if a >= 0 else -a
  else:
    return -a if a >= 0 else a


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

  #h = np.sign(h1) * (x2 - x1) # This is wrong. I made a worng conversion from C++ to python !
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

    if (x+h-x2)*(x+h-x1) > 0.0:
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




