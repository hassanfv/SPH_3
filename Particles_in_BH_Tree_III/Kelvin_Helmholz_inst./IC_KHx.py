
import numpy as np
import matplotlib.pyplot as plt
import random
import numpy as np

np.random.seed = 42

# Simulation parameters
domain_size = (1.0, 1.0, 1.0)  # Domain size (x, y, z)
resolution = 100  # Particles per unit length in one direction

# Calculate the number of particles based on the resolution and domain size
num_particles_x = int(resolution * domain_size[0])
num_particles_y = int(resolution * domain_size[1])
num_particles_z = int(resolution * domain_size[2])

# Initialize particle positions
x, y, z = np.meshgrid(
    np.linspace(0, domain_size[0], num_particles_x),
    np.linspace(0, domain_size[1], num_particles_y),
    np.linspace(0, domain_size[2], num_particles_z)
)

# Flatten the arrays
x = x.flatten()
y = y.flatten()
z = z.flatten()

# Initialize velocities and densities
mass = np.zeros_like(x)

N = Npart = x.shape[0]
Volume = 1.0 * 1.0 * 1.0

h = (Volume / Npart)**(1./3.)
h = h + np.zeros_like(x)

# Adding randomness to the positions
x = np.array([(tmp - 0.5 + h[0] * (-1.0 + random.random())) for tmp in x])
y = np.array([(tmp - 0.5 + h[0] * (-1.0 + random.random())) for tmp in y])
z = np.array([(tmp - 0.5 + h[0] * (-1.0 + random.random())) for tmp in z])


V_Low = 0.25
m_tot_L = 1.0 * V_Low
msph_L = m_tot_L / N / 4

V_High = 0.5
m_tot_H = 2.0 * V_High
msph_H = m_tot_H / N / 2

mass = mass + msph_L
mass[np.abs(y - 0.0) <= 0.25] = msph_H

#plt.figure(figsize=(13, 10))
#scatter = plt.scatter(x, y, s = 0.1, c = np.log10(mass), cmap='rainbow')
#plt.colorbar(scatter, label='np.log10(rho)')
#plt.show()

print(mass)


# Generate Initial Conditions - opposite moving streams with perturbation
w0 = 0.1
sigma = 0.05/np.sqrt(2.)
rho = 1. + (np.abs(y-0.0) < 0.25)
vx = -0.5 + (np.abs(y-0.0)<0.25)
vy = w0*np.sin(4*np.pi*x) * ( np.exp(-(y-0.25)**2/(2 * sigma**2)) + np.exp(-(y-0.75)**2/(2*sigma**2)) )

vz = 0.0 * vx

P = 2.5 * np.ones(N)

gamma = 5.0 / 3.0
u = P / (gamma - 1.0) / rho

print('u = ', u)


epsilon = h

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

mass = mass.astype(np.float32)
u = u.astype(np.float32)
h = h.astype(np.float32)
epsilon = epsilon.astype(np.float32)

# Save the arrays to a binary file:
N_tot = Npart
num = str(int(np.floor(N_tot/1000)))
filename = 'IC_KHX_' + num + 'k.bin'
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

print()
print("mass = ", mass)

plt.scatter(x, z, s = 0.5)
plt.show()



