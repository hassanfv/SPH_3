
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation #to make animation 

def lorenz(t, xyz, sigma, rho, beta):

  x, y, z = xyz

  dxyz = [sigma * (y - x),
          x * (rho - z) - y,
          x*y - beta*z
         ]

  return np.array(dxyz)



def rk4(t, xyz, dt, sigma, rho, beta):

  f1 = lorenz(t, xyz, sigma, rho, beta)

  f2 = lorenz(t + 0.5*dt, xyz + 0.5*dt*f1, sigma, rho, beta)

  f3 = lorenz(t + 0.5*dt, xyz + 0.5*dt*f2, sigma, rho, beta)

  f4 = lorenz(t + dt, xyz + dt*f3, sigma, rho, beta)

  xyz += (dt*f1 + 2*dt*f2 + 2*dt*f3 + dt*f4) / 6.0

  return xyz


#initial conditions
xyz0 = [1, 1, 1]
sigma, rho, beta = 10, 28, 8/3
t, dt = 0, 0.01
#integration loop
x_list = []
y_list = []
z_list = []
for i in range(10000):
    xyz = rk4(t, xyz0, dt, sigma, rho, beta)
    x_list.append(xyz[0])
    y_list.append(xyz[1])
    z_list.append(xyz[2])
    xyz0 = xyz
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




