
import numpy as np
import matplotlib.pyplot as plt


def traject(t):
  
  x = 1.5 * t
  return x**(1./3.) + np.sqrt(t) + 1.0

def deriv(t):
    return (1/3)*(1.5)**(1/3) * t**(-2/3) + 0.5 * t**(-0.5)



dt = 0.5

def FEuler(x_k, t_k):
  x_k_plus_1 = x_k + dt * deriv(t_k)
  return x_k_plus_1


t_k = 1.0
x_k = traject(t_k) # initial value

res = []
res.append([t_k, x_k])


while t_k <= 10:
  x_k = FEuler(x_k, t_k)
  t_k = t_k + dt

  res.append([t_k, x_k])

  print(f"At t = {t_k} the particle is at x = {x_k}")

res = np.array(res)

t = np.arange(1., 10, 0.01)
plt.plot(t, traject(t), color = 'k')
plt.plot(res[:, 0], res[:, 1], color = 'r', linestyle = "--")
plt.show()



