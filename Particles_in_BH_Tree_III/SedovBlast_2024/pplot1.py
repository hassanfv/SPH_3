
import numpy as np
import matplotlib.pyplot as plt
import struct

filename = 'KH-0.320000.bin'

def readArraysFromBinary(filename):
    with open(filename, 'rb') as file:
        # Read N
        N = np.fromfile(file, dtype=np.int32, count=1)[0]

        # Read the arrays from the file
        Typ = np.fromfile(file, dtype=np.int32, count=N)
        x = np.fromfile(file, dtype=np.float32, count=N)
        y = np.fromfile(file, dtype=np.float32, count=N)
        z = np.fromfile(file, dtype=np.float32, count=N)
        vx = np.fromfile(file, dtype=np.float32, count=N)
        vy = np.fromfile(file, dtype=np.float32, count=N)
        vz = np.fromfile(file, dtype=np.float32, count=N)
        rho = np.fromfile(file, dtype=np.float32, count=N)
        h = np.fromfile(file, dtype=np.float32, count=N)
        u = np.fromfile(file, dtype=np.float32, count=N)
        mass = np.fromfile(file, dtype=np.float32, count=N)

    return x, y, z, vx, vy, vz, rho, h, u, mass, Typ

# Usage
x, y, z, vx, vy, vz, rho, h, u, mass, Typ = readArraysFromBinary(filename)

print('Typ == 0 ===> ', np.sum(Typ == 0))
print()

n = np.where(u != 0.0)[0]
rho = rho[n]
u = u[n]
mass = mass[n]

print('sort(mass) = ', np.sort(mass))

h = h[n]
print('sort(h) = ', np.sort(h))
print()
print('sort(rho) = ', np.sort(rho))
print()

x = x[n]
y = y[n]
z = z[n]

rr = np.sqrt(x*x + y*y + z*z)


vx = vx[n]
vy = vy[n]
vz = vz[n]

#nz = np.where((rr < 0.50) & (np.abs(z) < 0.005))[0]
nz = np.where(np.abs(z) < 0.20)[0]

x = x[nz]
y = y[nz]
z = z[nz]

u = u[nz]
rho = rho[nz]

rr = (x*x + y*y + z*z)**0.5

h = h[nz]

vx = vx[nz]
vy = vy[nz]
vz = vz[nz]

vv = np.sqrt(vx*vx + vy*vy + vz*vz)

print('HHHH')

plt.figure(figsize=(13, 10))

scatter = plt.scatter(x, y, s = 1.0, c = np.log10(rho), cmap='rainbow')

plt.colorbar(scatter, label='np.log10(rho)')

xy = 5.1

plt.xlim(-xy, xy)
plt.ylim(-xy, xy)
#plt.xlabel('X')
#plt.ylabel('Y')

# Optional: Adjust layout
plt.tight_layout()


plt.savefig('XSB.png')

plt.show()




