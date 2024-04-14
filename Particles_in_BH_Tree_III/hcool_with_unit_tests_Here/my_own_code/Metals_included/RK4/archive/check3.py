
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation #to make animation 


sigma = 10.0
beta = 8.0 / 3.0
rho = 28.0

print('check!')

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
y0 = [-8.0, 8.0, 27.0]

t, dt = 0, 0.01
#integration loop
x_list = []
y_list = []
z_list = []
for i in range(10000):
    xyz = rk4SingleStep(lorenz, dt, t, y0)
    x_list.append(xyz[0])
    y_list.append(xyz[1])
    z_list.append(xyz[2])
    y0 = xyz
    t += dt


#plotting
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x_list, y_list, z_list, lw=.5)
ax.set_title("Lorenz System Trajectory")
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
plt.show()




