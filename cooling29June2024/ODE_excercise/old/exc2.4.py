
import numpy as np
import matplotlib.pyplot as plt


t_n = 0.
u_n = 1.
v_n = 0.

h = 0.1

t_n_plus_1 = t_n + h

for i in range(2):
  
  u_n_plus_1 = u_n + h * (-2. * u_n + v_n)
  v_n_plus_1 = v_n + h * (-u_n - 2. * v_n)

  print(t_n_plus_1, u_n_plus_1)
  print(t_n_plus_1, v_n_plus_1)

  t_n_plus_1 += h
  u_n = u_n_plus_1
  v_n = v_n_plus_1

