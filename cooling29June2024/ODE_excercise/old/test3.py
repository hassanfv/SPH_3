
import numpy as np
import matplotlib.pyplot as plt

#---- func
def func(t_n, x_n):
  return (1. - 2. * t_n) * x_n


t_n = 0.
x_n = 1.

h = 0.1

t = np.arange(0, 4., 0.001)
plt.plot(t, np.exp((1./4. - (t - 1./2.)**2)), color = 'orange')


#------ TS(1) - Taylor Series order 1 -----
t_res = []
x_res = []
t_res.append(t_n)
x_res.append(x_n)
while t_n <= 4.0:
  
  x_n_p_1 = x_n + h * func(t_n, x_n)

  t_n += h
  t_res.append(t_n)
  x_res.append(x_n_p_1)
  
  x_n = x_n_p_1
#-------------------------------------------
plt.scatter(t_res, x_res, s = 20, color = 'k', label = 'TS(1)')


#------ TS(1) - Taylor Series order 1 -----
t_n = 0.
x_n = 1.
t_res = []
x_res = []
t_res.append(t_n)
x_res.append(x_n)
while t_n <= 4.0:
  
  x_n_p_1 = x_n + h * func(t_n, x_n)
  t_n_p_1 = t_n + h
  
  x_n_p_1 = x_n + 0.5 * h * (func(t_n_p_1, x_n_p_1) + func(t_n, x_n))

  t_res.append(t_n_p_1)
  x_res.append(x_n_p_1)
  
  t_n = t_n_p_1
  x_n = x_n_p_1
#-------------------------------------------
plt.scatter(t_res, x_res, s = 20, color = 'lime', label = 'Trapez.')

print(x_res)



plt.legend()
plt.show()
  
