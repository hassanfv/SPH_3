

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



