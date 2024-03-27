
import numpy as np
import matplotlib.pyplot as plt
import random

#===== createGrid
def createGrid(size):

  step = 1.0 / (Grid / 2.0)

  x = np.zeros(size)
  y = np.zeros(size)
  z = np.zeros(size)

  for i in range(size):
    ix = i / (Grid * Grid)
    iy = (i / Grid) % Grid
    iz = i % Grid
  
    x[i] = (ix - (Grid / 2 - 0.5)) * step
    y[i] = (iy - (Grid / 2 - 0.5)) * step
    z[i] = (iz - (Grid / 2 - 0.5)) * step
  
  return x, y, z
  

Grid = 161; # size of the grid
Mtot = 1.0; # total mass of the sphere

size = Grid * Grid * Grid

x, y, z = createGrid(size)
  
# Stretching initial conditions to get 1/r density distribution
for i in range(size):
  r = np.sqrt(np.sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]))
  if (r > 0):
    x[i] *= r
    y[i] *= r
    z[i] *= r




# Counting particles inside the unit sphere
Npart = 0
for i in range(size):
  if (np.sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]) < 1.0):
    Npart += 1

xx = np.zeros(Npart)
yy = np.zeros(Npart)
zz = np.zeros(Npart)

vx = np.zeros(Npart)
vy = np.zeros(Npart)
vz = np.zeros(Npart)

eps = np.zeros(Npart)

Typ = np.zeros(Npart)

mass = np.zeros(Npart)
Uthermal = np.zeros(Npart)

particle_mass = Mtot / Npart


pert = 0.005

k = 0;
for i in range(size):
  if (np.sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]) < 1.0):
      xx[k] = x[i] + pert * (-1.0 + 2.0 * random.random())
      yy[k] = y[i] + pert * (-1.0 + 2.0 * random.random())
      zz[k] = z[i] + pert * (-1.0 + 2.0 * random.random())
      
      vx[k] = 0.0
      vy[k] = 0.0
      vz[k] = 0.0
      
      Typ[k] = 0
      
      mass[k] = particle_mass
      Uthermal[k] = 0.05
      k += 1

  
print(f"We use {Npart} particles.")

R = 1.0  # radius of the sphere
h = np.ones(Npart) * 2.0 * R / Npart ** (1.0 / 3.0)

eps = h.copy()

u = Uthermal.copy()

#------------------------------------------------------------------------
x = xx
y = yy
z = zz
print(x.shape)

#vx = v[:, 0]
#vy = v[:, 1]
#vz = v[:, 2]

epsilon = h.copy()

Typ = np.zeros_like(x)

x = np.round(x, 5)
y = np.round(y, 5)
z = np.round(z, 5)

#vx = np.round(vx, 5)
#vy = np.round(vy, 5)
#vz = np.round(vz, 5)

Typ = Typ.astype(np.int32)

x = x.astype(np.float32)
y = y.astype(np.float32)
z = z.astype(np.float32)

vx = vx.astype(np.float32)
vy = vy.astype(np.float32)
vz = vz.astype(np.float32)

mass = mass.astype(np.float32)
u = u.astype(np.float32)
h = h.astype(np.float32)
epsilon = epsilon.astype(np.float32)

print('sort h = ', np.sort(h))

# Save the arrays to a binary file:
N_tot = Npart = len(x)
num = str(int(np.floor(N_tot/1000)))
filename = 'IC_Evrard24_' + num + 'k.bin'
with open(filename, "wb") as file:
  # Save the rest of the arrays:
  file.write(Typ.tobytes())
  file.write(x.tobytes())
  file.write(y.tobytes())
  file.write(z.tobytes())

  file.write(vx.tobytes())
  file.write(vy.tobytes())
  file.write(vz.tobytes())

  file.write(mass.tobytes())
  file.write(h.tobytes())
  file.write(epsilon.tobytes())
  file.write(u.tobytes())


plt.scatter(x, y, s = 0.1, color = 'black')

plt.savefig('IC.png')

plt.show()


