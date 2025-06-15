
# RK4 on Lorenz 1963 - 3D plot !

import numpy as np
import matplotlib.pyplot as plt
import time


#----- rk4singlestep
def rk4singlestep(fun, dt, t0, y0):

  f1 = fun(t0, y0)
  f2 = fun(t0 + dt / 2., y0 + (dt / 2.) * f1)
  f3 = fun(t0 + dt / 2., y0 + (dt / 2.) * f2)
  f4 = fun(t0 + dt, y0 + dt * f3)
  
  yout = y0 + (dt / 6.) * (f1 + 2. * f2 + 2. * f3 + f4)
  
  return yout


#----- lorenz
def lorenz(t, y):
  
  x_dot = sigma * (y[1] - y[0])
  y_dot = y[0] * (rho - y[2]) - y[1]
  z_dot = y[0] * y[1] - beta * y[2]
  
  return np.array([x_dot, y_dot, z_dot])


sigma = 10
beta = 8. / 3.
rho = 28.

# Initial condition
y0 = [-8., 8, 27.]

dt = 0.01

tin = 0
yin = y0

traj = [yin]
t = [tin]

TA = time.time()

while tin <= 10:

  yin = rk4singlestep(lorenz, dt, tin, yin)
  
  tin += dt
  
  traj.append(list(yin))
  t.append(tin)

print("Elapsed time in main loop = ", time.time() - TA)

traj = np.array(traj)

print(traj)
print()
print(traj.shape)

#plt.scatter(traj[:, 1], traj[:, 2], s = 3, color = 'k')
#plt.show()

ax = plt.figure().add_subplot(projection = '3d')
ax.scatter(traj[:, 0], traj[:, 1], traj[:, 2], s = 3, color = 'k')
ax.plot(traj[:, 0], traj[:, 1], traj[:, 2], color = 'gray')

plt.show()


