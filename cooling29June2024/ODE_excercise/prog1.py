
import numpy as np
import matplotlib.pyplot as plt


def func(t_n, x_n):
  return (1. - 2. * t_n) * x_n 


t_n = 0.
x_n = 1.

h = 0.2

res = []

res.append([t_n, x_n])

while t_n < 4.:
  
  x_n_p_1 = x_n + h * func(t_n, x_n)
  
  t_n_p_1 = t_n + h
  
  res.append([t_n_p_1, x_n_p_1])
  
  t_n = t_n_p_1
  x_n = x_n_p_1


res = np.array(res)

tt = res[:, 0]
xx = res[:, 1]  

t = np.arange(0, 4, 0.001)
x = np.exp(1./4. - (t - 1./2.)**2)

plt.plot(t, x, color = 'k')
plt.scatter(tt, xx, color = 'orange')


plt.show()

