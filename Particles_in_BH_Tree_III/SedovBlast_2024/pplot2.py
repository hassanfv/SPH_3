
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

vx = vx[n]
vy = vy[n]
vz = vz[n]

#nz = np.where((rr < 0.50) & (np.abs(z) < 0.005))[0]
nz = np.where(np.abs(z) < 1110.20)[0]

x = x[nz]
y = y[nz]
z = z[nz]

u = u[nz]
rho = rho[nz]

h = h[nz]

vx = vx[nz]
vy = vy[nz]
vz = vz[nz]

vv = np.sqrt(vx*vx + vy*vy + vz*vz)

rr = np.sqrt(x*x + y*y + z*z)

v_radial = (vx * x + vy * y + vz * z) / rr

rgrid = np.linspace(0, 3.0, 100)

print(rgrid)


rho_bin = np.zeros_like(rgrid)
vel_bin = np.zeros_like(rgrid)

res = []

for i in range(len(rgrid) - 1):
  
  nt = np.where((rr >= rgrid[i]) & (rr < rgrid[i+1]))[0]
  print(nt)
  
  rho_bin[i] = np.mean(rho[nt])
  vel_bin[i] = np.mean(v_radial[nt])
  
  res.append([0.5*(rgrid[i]+rgrid[i+1]), rho_bin[i], vel_bin[i]])


res = np.array(res)
rx = res[:, 0]
rhox = res[:, 1]
vel = res[:, 2]


plt.figure(figsize=(15, 7))

#------> Model <--------
rho_data = np.loadtxt("rho.dat")
r_s_rho, rho_s = rho_data[:, 0], rho_data[:, 1]

v_data = np.loadtxt("v.dat")
r_s_v, v_s = v_data[:, 0], v_data[:, 1]
#-----------------------


plt.subplot(1, 2, 1)
plt.scatter(rr, rho, s = 0.05, color = 'orange')
plt.scatter(rx, rhox, s = 20, color = 'b')
plt.plot(r_s_rho, rho_s, color = 'k')
plt.xlim(0, 3)


plt.subplot(1, 2, 2)
plt.scatter(rr, v_radial, s = 0.05, color = 'orange')
plt.scatter(rx, vel, s = 20.0, color = 'b')
plt.plot(r_s_v, v_s, color = 'k')
plt.xlim(0, 3)


plt.savefig('rho_vs_r.png')

plt.show()




