import numpy as np
import matplotlib.pyplot as plt
import time

#----- rk4
def rk4(y, dydx, x, h):

  n = len(y)
  
  yout = np.zeros(n)
  dyt = np.zeros(n)
  yt = np.zeros(n)
  dym = np.zeros(n)
  
  hh = 0.5 * h
  h6 = h / 6.0
  
  xh = x + hh # i.e. x + h/2
  
  dydx = derivs(x, y) # this is k1
  
  for i in range(n):
    yt[i] = y[i] + hh * dydx[i] # y0 + h/2 * k1
  dyt = derivs(xh, yt) # this is k2.
  
  for i in range(n):
    yt[i] = y[i] + hh * dyt[i]
  dym = derivs(xh, yt) # this is k3.
  
  for i in range(n):
    yt[i] = y[i] + h * dym[i]
    dym[i] += dyt[i] # clever idea !
  dyt = derivs(x + h, yt)
  
  for i in range(n):
    yout[i] = y[i] + h6 * (dydx[i] + dyt[i] + 2.0 * dym[i])
  
  return yout
  


def derivs(x, y):
    sigma = 10.0
    rho = 28.0
    beta = 8.0 / 3.0

    dydx = np.zeros_like(y)
    dydx[0] = sigma * (y[1] - y[0])              # dx/dt = sigma*(y - x)
    dydx[1] = y[0] * (rho - y[2]) - y[1]         # dy/dt = x*(rho - z) - y
    dydx[2] = y[0] * y[1] - beta * y[2]          # dz/dt = x*y - beta*z
    return dydx

h = 0.01  # Step size
x = 0.0  # Initial time
y = np.array([1.0, 1.0, 1.0])  # Initial values for x, y, z

dydx = derivs(x, y)

res = []

TA = time.time()
for j in range(10000):  # Running the RK4 for more steps to observe the chaotic behavior
    y = rk4(y, dydx, x, h)
    res.append([x, y[0], y[1], y[2]])
    dydx = derivs(x, y)
    x += h

print('Elapsed time = ', time.time() - TA)

res = np.array(res)

# Plotting the results
plt.figure(figsize=(12, 8))
plt.plot(res[:, 0], res[:, 1], label='X(t)')
plt.plot(res[:, 0], res[:, 2], label='Y(t)')
plt.plot(res[:, 0], res[:, 3], label='Z(t)')
plt.legend()
plt.title('Lorenz Attractor')
plt.xlabel('Time')
plt.ylabel('States')
plt.grid(True)
plt.show()




