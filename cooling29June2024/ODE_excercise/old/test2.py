
import numpy as np
import matplotlib.pyplot as plt

t_n = 0.
u_n = 2.
v_n = -9.

h = 0.0001

t_res = []
x_res = []

t_res.append(t_n)
x_res.append(u_n)

for i in range(10000):

  u_n_p_1 = u_n + h * v_n
  v_n_p_1 = v_n + h * 50 * np.exp(-5. * t_n)
  
  t_res.append(t_n)
  x_res.append(u_n_p_1) # Note: u = x
  
  u_n = u_n_p_1
  v_n = v_n_p_1
  
  t_n += h

t = np.arange(0., 1, 0.001)

plt.plot(t, t+2.*np.exp(-5. * t), color = 'r')
plt.scatter(t_res, x_res, s = 10, color = 'k')


plt.show()


