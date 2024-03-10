
import numpy as np
import matplotlib.pyplot as plt
import struct

filename = './Outputs/G-0.011400.bin'


tCode = np.float64(filename[-12:-4]) / 10.0 #!!!!!!!!!!!!!!!!! Double check if during the file-saving the multiplication by 10 is actually done!!!!!!!!!!


unit_velocity_cgs = 3.1338e+06 # cm/s #!!!!!!!!!!!!!!!!!!!!!!!!
unit_u = 9.82068e+12 #!!!!!!!!!!!!!!!!!!!!!!!!
unit_rho = 1.54506e-23 # !!!!!!!!!!!!!!!!!!!
unitTime_in_s = 9.84748e+14 #!!!!!!!!!!!!!!!!!!!!!!!!!

t_in_kyrs = tCode * unitTime_in_s / 3600./24./365.25/1000.

def readBinaryFile(filename):
    with open(filename, 'rb') as f:
        # Read N and N_ionFrac
        N, N_ionFrac = struct.unpack('ii', f.read(2 * 4))  # 4 bytes each for two integers

        # Read arrays
        Typ = np.array(struct.unpack(f'{N}i', f.read(N * 4)))
        x = np.array(struct.unpack(f'{N}f', f.read(N * 4)))
        y = np.array(struct.unpack(f'{N}f', f.read(N * 4)))
        z = np.array(struct.unpack(f'{N}f', f.read(N * 4)))
        vx = np.array(struct.unpack(f'{N}f', f.read(N * 4)))
        vy = np.array(struct.unpack(f'{N}f', f.read(N * 4)))
        vz = np.array(struct.unpack(f'{N}f', f.read(N * 4)))
        rho = np.array(struct.unpack(f'{N}f', f.read(N * 4)))
        h = np.array(struct.unpack(f'{N}f', f.read(N * 4)))
        u = np.array(struct.unpack(f'{N}f', f.read(N * 4)))
        mass = np.array(struct.unpack(f'{N}f', f.read(N * 4)))
        ionFrac = np.array(struct.unpack(f'{N_ionFrac}f', f.read(N_ionFrac * 4)))

    # Return the data
    return N, N_ionFrac, Typ, x, y, z, vx, vy, vz, rho, h, u, mass, ionFrac 

# Usage
N, N_ionFrac, Typ, x, y, z, vx, vy, vz, rho, h, u, mass, ionFrac = readBinaryFile(filename)

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

nz = np.where((rr < 0.50) & (np.abs(z) < 0.005))[0]
#nz = np.where(np.abs(z) < 0.03)[0]

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

kB = 1.3807e-16
mu = 0.61
mH = 1.673534e-24
gamma = 5./3.
Temp = (gamma - 1) * mH / kB * mu * u * unit_u
print('sort T = ', np.sort(Temp))#[-5:])

rho_cgs = rho * unit_rho
XH = 0.7
nH_cgs = rho_cgs * XH / mH


unit_density_in_cgs = unit_rho

nH = rho * unit_density_in_cgs * XH /mH

print(f'max(nH) = {max(rho * unit_density_in_cgs * XH /mH)} ---->  note this could be different from the global Max(nH).')
print()

print()
print('sort(Temp) = ', np.sort(Temp))
print()

nk = np.where((Temp > 1e6) & (Temp < 1e9))[0]
xT = x[nk]
yT = y[nk]
zT = z[nk]
TempT = Temp[nk]
#print('nk = ', nk)
print()
print(f'Number of particles with 1e6 < Temp < 1e9 is {len(TempT)}')
#print('Temp[nk] = ', Temp[nk])
print()

ntmp = np.where((Temp > 60000) & (Temp < 100000))[0]
#print('ntmp = ', ntmp)

print()
print('Number of particles with T > 1e6 K (Note that this is for the thing layer) = ', len(np.where(Temp > 1e6)[0]))
print()
print('Number of particles with T > 1e9 K (Note that this is for the thing layer) = ', len(np.where(Temp > 1e9)[0]))
print()

plt.figure(figsize=(13, 10))

#scatter = plt.scatter(x, y, c=np.log10(Temp), cmap='rainbow', s=0.2)
scatter = plt.scatter(x, y, c=np.log10(nH_cgs), cmap='rainbow', s=0.2)

plt.colorbar(scatter, label='log10(Temperature)')

xy = 0.55

plt.xlim(-xy, xy)
plt.ylim(-xy, xy)
plt.xlabel('X')
plt.ylabel('Y')

# Optional: Adjust layout
plt.tight_layout()


plt.savefig('RT.png')

plt.show()




