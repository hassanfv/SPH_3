
import numpy as np
import matplotlib.pyplot as plt
import random
import numpy as np

np.random.seed = 42

# Simulation parameters
domain_size = (1.0, 1.0, 1.0)  # Domain size (x, y, z)
resolution = 100  # Particles per unit length in one direction
perturbation_amplitude = 0.02 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
velocity_shear = 0.5 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
mass_1 = 1.0  # mass for the first fluid #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
mass_2 = 2.0  # mass for the second fluid (denser) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
vx = np.zeros_like(x)
vy = np.zeros_like(y)
vz = np.zeros_like(z)
mass = np.zeros_like(x)

Npart = x.shape[0]
Volume = 1.0 * 1.0 * 1.0

h = (Volume / Npart)**(1./3.)
h = h + np.zeros_like(x)

# Adding randomness to the positions
x = np.array([(tmp - 0.5 + h[0] * (-1.0 + random.random())) for tmp in x])
y = np.array([(tmp - 0.5 + h[0] * (-1.0 + random.random())) for tmp in y])
z = np.array([(tmp - 0.5 + h[0] * (-1.0 + random.random())) for tmp in z])


u0 = 0.1 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
u = np.zeros_like(x) + u0

epsilon0 = 0.005 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
epsilon = np.zeros_like(x) + epsilon0

Typ = np.zeros_like(x)

# Apply velocity shear, mass contrast, and perturbation
for i in range(len(x)):
    if y[i] > domain_size[1] / 2:
        vx[i] = velocity_shear
        mass[i] = mass_1
    else:
        vx[i] = -velocity_shear
        mass[i] = mass_2
    # Add perturbation at the interface
    y[i] += perturbation_amplitude * np.sin(6 * np.pi * x[i] / domain_size[0])


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


plt.scatter(x, y, s = 0.5)
plt.show()



