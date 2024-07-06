
import numpy as np


def f(t_n, x_n):
  return 2. * x_n * (1. - x_n)

t_n = 10.
h = 0.2
x_n = 1. / 5.

t_n += h

while (t_n < 11.):
  x_nPlus1 = x_n + h * f(t_n, x_n)
  print(t_n, x_nPlus1)
  
  x_n = x_nPlus1
  t_n += h



