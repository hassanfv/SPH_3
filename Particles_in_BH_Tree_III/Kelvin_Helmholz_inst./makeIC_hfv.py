################################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
#               2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
################################################################################

import h5py
import numpy as np
import matplotlib.pyplot as plt
import random

# Generates a swift IC file for the Kelvin-Helmholtz vortex in a periodic box

# Parameters
gamma = 5.0 / 3.0  # Gas adiabatic index
P0 = 2.5  # Pressure
rho0 = 1.0  # Density
d = 0.0317  # Thickness of the transition layer
B = 0.0005  # Amplitude of the seed velocity
N1D = 96  # Number of particles in one dimension
# ---------------------------------------------------

N = N1D ** 3
x = np.linspace(0.0, 1.0, N1D + 1)
x = 0.5 * (x[1:] + x[:-1])
y = x
z = x
xx, yy, zz = np.meshgrid(x, y, z)
pos = np.zeros((N, 3))
pos[:, 0] = xx.reshape((N))
pos[:, 1] = yy.reshape((N))
pos[:, 2] = zz.reshape((N))

L = N1D
boxSize = 1.0
pert = 0.1

for i in range(L):
    for j in range(L):
        for k in range(L):
            index = i*L*L + j*L + k
            x = i * boxSize / L + boxSize / (2*L)
            y = j * boxSize / L + boxSize / (2*L)
            z = k * boxSize / L + boxSize / (2*L)

            pos[index,0] = x - 0.5 + random.random() * pert * boxSize/(2.*L)
            pos[index,1] = y - 0.5 + random.random() * pert * boxSize/(2.*L)
            pos[index,2] = z - 0.5 + random.random() * pert * boxSize/(2.*L)
#--------------------------------------------------



h = np.ones(N) * 2.0 / N1D

vol = 1.0

# Generate extra arrays
v = np.zeros((N, 3))
ids = np.linspace(1, N, N)
m = np.ones(N) * rho0 * vol / N
u = np.ones(N) * P0 / (rho0 * (gamma - 1.0))

v[pos[:, 1] <= 0.25 - d, 0] = -0.5
v[(pos[:, 1] < 0.25 + d) & (pos[:, 1] > 0.25 - d), 0] = (
    -0.5
    + 0.5 * (pos[(pos[:, 1] < 0.25 + d) & (pos[:, 1] > 0.25 - d), 1] + d - 0.25) / d
)
v[(pos[:, 1] <= 0.75 - d) & (pos[:, 1] >= 0.25 + d), 0] = 0.5
v[(pos[:, 1] < 0.75 + d) & (pos[:, 1] > 0.75 - d), 0] = (
    0.5 - 0.5 * (pos[(pos[:, 1] < 0.75 + d) & (pos[:, 1] > 0.75 - d), 1] + d - 0.75) / d
)
v[pos[:, 1] >= 0.75 + d, 0] = -0.5

v[:, 1] = (
    B
    * np.sin(4.0 * np.pi * pos[:, 0])
    * (
        np.exp(-(pos[:, 1] - 0.25) ** 2 / 32.0 / d ** 2)
        + np.exp(-(pos[:, 1] - 0.75) ** 2 / 32.0 / d ** 2)
    )
)

#------------------------------------------------------------------------

x = pos[:, 0]
y = pos[:, 1]
z = pos[:, 2]
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
filename = 'IC_KH_' + num + 'k.bin'
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





