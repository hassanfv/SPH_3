
# Implicit

import numpy as np
import matplotlib.pyplot as plt


def func(t_n, x_n):
  return (1. - 2. * t_n) * x_n 


tn = 0.
xn = 1.

h = 0.2

tn +=h

res = []

res.append([tn, xn])

xkp1 = xn
xk = xkp1

for i in range(10):

  g_xnp1 = xk - xn - h * func(tn, xk)
  deriv_g_xnp1 = 1 - h * (1. - 2. * tn)

  xkp1 = xk - g_xnp1 / deriv_g_xnp1

  ratio = np.abs(xk - xkp1) / xk

  xk = xkp1

  print(xk, ratio, g_xnp1, deriv_g_xnp1, xkp1)

