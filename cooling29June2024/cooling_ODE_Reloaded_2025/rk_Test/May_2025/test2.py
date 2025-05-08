
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

def RK2(x_k, t_k):
  x_k_plus_1 = x_k + dt * deriv(t_k + dt/2)
  return x_k_plus_1


#--- FEuler ----
t_k = 1.0
x_k = traject(t_k) # initial value
resFE = []
resFE.append([t_k, x_k])
while t_k <= 10:
  x_k = FEuler(x_k, t_k)
  t_k = t_k + dt
  resFE.append([t_k, x_k])

#---- RK2 ----
t_k = 1.0
x_k = traject(t_k) # initial value
resRK2 = []
resRK2.append([t_k, x_k])
while t_k <= 10:
  x_k = RK2(x_k, t_k)
  t_k = t_k + dt
  resRK2.append([t_k, x_k])

resFE = np.array(resFE)
resRK2 = np.array(resRK2)

t = np.arange(1., 10, 0.01)
plt.plot(t, traject(t), color = 'k', label = 'Actual trajectory', linewidth = 5)
plt.plot(resFE[:, 0], resFE[:, 1], color = 'r', linestyle = "--", label = 'FE')
plt.plot(resRK2[:, 0], resRK2[:, 1], color = 'lime', linestyle = "--", label = 'RK2')

plt.legend()

plt.show()



