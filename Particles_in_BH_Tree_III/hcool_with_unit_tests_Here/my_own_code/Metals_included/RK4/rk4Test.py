
#Ref: https://www.youtube.com/watch?v=vNoFdtcPFdk

import numpy as np
import matplotlib.pyplot as plt


sigma = 10.0
beta = 8.0 / 3.0
rho = 28.0


#----- lorenz
def lorenz(t, y):

  dy = [
         sigma * (y[1] - y[0]),
         y[0] * (rho - y[2]) - y[1],
         y[0] * y[1] - beta * y[2]
       ]
  
  return np.array(dy)



#===== rk4SingleStep
def rk4SingleStep(func, dt, t0, y0):

  f1 = func(t0, y0)
  f2 = func(t0 + dt/2.0, y0 + (dt/2.0) * f1)
  f3 = func(t0 + dt/2.0, y0 + (dt/2.0) * f2)
  f4 = func(t0 + dt, y0 + dt * f3)
  
  yout = y0 + (dt / 6.0) * (f1 + 2.0 * f2 + 2.0 * f3 + f4)
  
  return yout


# Initial condition
y0 = np.array([-8.0, 8.0, 27.0])

# Compute trajectory

dt = 0.01
T = 100.0

num_time_pts = int(T / dt)
t = np.linspace(0.0, T, num_time_pts)


Y = np.zeros((3, num_time_pts))
Y[:, 0] = y0
yin = y0

for i in range(num_time_pts - 1):

  yout = rk4SingleStep(lorenz, dt, t[i], yin)
  Y[:, i+1] = yout
  yin = yout


print(Y)


ax = plt.figure().add_subplot(projection = '3d')
ax.plot(Y[0, :], Y[1, :], Y[2, :], 'b', lw = 0.5)
plt.show()










