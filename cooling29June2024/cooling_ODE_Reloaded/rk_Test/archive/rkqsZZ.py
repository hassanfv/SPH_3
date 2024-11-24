

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

  # When we exit the WHILE loop, we have the "y" which is evolved for one step h with the best h (i.e. with a h that gives errmax <= 1.0)!
  if errmax > ERRCON:
    hnext = SAFETY * h * errmax ** PGROW
  else:
    hnext = 5.0 * h

  x += h
  hdid = h

  for i in range(n):
    y[i] = ytemp[i]

  return x, y, hdid, hnext




