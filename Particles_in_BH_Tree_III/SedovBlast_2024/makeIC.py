###############################################################################
 # This file is part of SWIFT.
 # Coypright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 #                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 # 
 # This program is free software: you can redistribute it and/or modify
 # it under the terms of the GNU Lesser General Public License as published
 # by the Free Software Foundation, either version 3 of the License, or
 # (at your option) any later version.
 # 
 # This program is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 # GNU General Public License for more details.
 # 
 # You should have received a copy of the GNU Lesser General Public License
 # along with this program.  If not, see <http://www.gnu.org/licenses/>.
 # 
 ##############################################################################

import random
import numpy as np
import matplotlib.pyplot as plt

# Generates a swift IC file for the Sedov blast test in a periodic cubic box

# Parameters
periodic= 1      # 1 For periodic box
boxSize = 10.
L = 131           # Number of particles along one axis
rho = 1.          # Density
P = 1.e-5         # Pressure
E0= 1.e2          # Energy of the explosion
pert = 0.1
gamma = 5./3.     # Gas adiabatic index

#---------------------------------------------------
numPart = L**3
mass = boxSize**3 * rho / numPart
internalEnergy = P / ((gamma - 1.)*rho)

if L%2 == 0:
    print("Number of particles along each dimension must be odd.")
    exit()

#Generate particles
coords = np.zeros((numPart, 3))
v      = np.zeros((numPart, 3))
m      = np.zeros((numPart, 1))
h      = np.zeros((numPart, 1))
u      = np.zeros((numPart, 1))
ids    = np.zeros((numPart, 1), dtype='L')

for i in range(L):
    for j in range(L):
        for k in range(L):
            index = i*L*L + j*L + k
            x = i * boxSize / L + boxSize / (2*L)
            y = j * boxSize / L + boxSize / (2*L)
            z = k * boxSize / L + boxSize / (2*L)
            coords[index,0] = x
            coords[index,1] = y
            coords[index,2] = z
            v[index,0] = 0.
            v[index,1] = 0.
            v[index,2] = 0.
            m[index] = mass
            h[index] = 2.251 / 2 * boxSize / L
            u[index] = internalEnergy
            ids[index] = index
            if np.sqrt((x - boxSize/2.)**2 + (y - boxSize/2.)**2 + (z - boxSize/2.)**2) < 2.01 * boxSize/L:
                u[index] = u[index] + E0 / (33. * mass)
                print("Particle " , index , " set to detonate.")
            coords[index,0] = x - 5.0 + random.random() * pert * boxSize/(2.*L)
            coords[index,1] = y - 5.0 + random.random() * pert * boxSize/(2.*L)
            coords[index,2] = z - 5.0 + random.random() * pert * boxSize/(2.*L)
#--------------------------------------------------

x = coords[:, 0]
y = coords[:, 1]
z = coords[:, 2]
print(x.shape)

vx = v[:, 0]
vy = v[:, 1]
vz = v[:, 2]

epsilon = h.copy()

Typ = np.zeros_like(x)

x = np.round(x, 5)
y = np.round(y, 5)
z = np.round(z, 5)

vx = np.round(vx, 5)
vy = np.round(vy, 5)
vz = np.round(vz, 5)

Typ = Typ.astype(np.int32)

x = x.astype(np.float32)
y = y.astype(np.float32)
z = z.astype(np.float32)

vx = vx.astype(np.float32)
vy = vy.astype(np.float32)
vz = vz.astype(np.float32)

mass = m.astype(np.float32)
u = u.astype(np.float32)
h = h.astype(np.float32)
epsilon = epsilon.astype(np.float32)

# Save the arrays to a binary file:
N_tot = Npart = len(x)
num = str(int(np.floor(N_tot/1000)))
filename = 'IC_SedovBlast_' + num + 'k.bin'
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
plt.show()



