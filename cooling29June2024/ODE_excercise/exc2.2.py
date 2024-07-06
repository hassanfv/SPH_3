
import numpy as np
import matplotlib.pyplot as plt


h = 1. / 20.
lamb = -10.

t_n = 0.0
x_n = 1.0

t_euler = []
x_euler = []

t_euler.append(t_n)
x_euler.append(x_n)

for i in range(10):

  x_n_plus_1 = x_n + h * lamb * x_n
  
  t_euler.append(t_n)
  x_euler.append(x_n_plus_1)
  
  t_n += h
  x_n = x_n_plus_1

print(x_euler)

t = np.arange(0, 1.0, 0.01)

plt.plot(t, np.exp(lamb * t), color = 'red')
plt.scatter(t_euler, x_euler, s = 10, color = 'b')


plt.show()

