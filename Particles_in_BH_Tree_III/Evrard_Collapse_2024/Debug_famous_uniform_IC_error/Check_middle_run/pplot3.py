
import numpy as np
import matplotlib.pyplot as plt
import struct
import pickle
import csv

filename = 'KH-8.004113.bin'

gamma = 5./3.

#loads 1D PPM results
def load_ppm_result():
	gamma = 5./3.
	rost = 3./4./np.pi
	est = 1.054811e-1  / 1.05
	pst = rost*est
	vst = np.sqrt(est)
	rst = 1.257607
	time = 0

	radius = np.zeros(350)
	rho = np.zeros(350)
	vr = np.zeros(350)
	press = np.zeros(350)

	with open('./ppm_profile/ppm1oaf') as csvfile:
		readCSV = csv.reader(csvfile)
		line = 0
		for row in readCSV:
			line = line+1
			values = row[0].split()
			if(line == 1):
				time = values[1]
				continue
			if(line == 352):
				break

			radius[line -2] = float(values[1]) /rst*1e-11
			rho[line -2] = float(values[2]) /rost
			vr[line -2] = float(values[4]) /vst*1e-8
			press[line -2] = float(values[3])/pst*1e-16

	rho=rho*(3.0/(4*np.pi))
	press = press*(3.0/(4*np.pi))

	entropy = press / rho**gamma

	return radius, rho, vr, entropy, press


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

print('hh = ', np.sort(h))


vx = vx[nz]
vy = vy[nz]
vz = vz[nz]

vv = np.sqrt(vx*vx + vy*vy + vz*vz)

rr = np.sqrt(x*x + y*y + z*z)

v_radial = (vx * x + vy * y + vz * z) / rr

#rgrid = np.linspace(0, 1.5, 100)
rgrid = 10**(np.linspace(np.log10(1e-3), np.log10(1.5), 100))

print(rgrid)


rho_bin = np.zeros_like(rgrid)
vel_bin = np.zeros_like(rgrid)
uTh = np.zeros_like(rgrid)

res = []

for i in range(len(rgrid) - 1):
  
  nt = np.where((rr >= rgrid[i]) & (rr < rgrid[i+1]))[0]
  print(nt)
  
  rho_bin[i] = np.mean(rho[nt])
  vel_bin[i] = np.mean(v_radial[nt])
  uTh[i] = np.mean(u[nt])
  
  res.append([0.5*(rgrid[i]+rgrid[i+1]), rho_bin[i], vel_bin[i], uTh[i]])


res = np.array(res)
rx = res[:, 0]
rhox = res[:, 1]
vel = res[:, 2]
ux = res[:, 3]


plt.scatter(x, y, s = 0.01)
plt.show()

#----------------------------------------
radius_ppm, rho_ppm, vr_ppm, entropy, press = load_ppm_result()

u_ppm = entropy * rho_ppm**(gamma - 1.0) / (gamma - 1.0) # see load_ppm_result() function above!
#----------------------------------------


plt.figure(figsize=(12, 8))

plt.subplot(2, 3, 1)
plt.scatter(rr, v_radial, s = 0.05, color = 'orange')
plt.scatter(rx, vel, s = 20.0, color = 'b')
plt.plot(radius_ppm, vr_ppm, linestyle = '--', color = 'k')
plt.xscale('log')
plt.xlim(5e-3, 1.25)
plt.ylim(-2, 2)
plt.title('radial velocity')


plt.subplot(2, 3, 2)
plt.scatter(rr, rho, s = 0.05, color = 'orange')
plt.scatter(rx, rhox, s = 20, color = 'b')
plt.plot(radius_ppm, rho_ppm, linestyle = '--', color = 'k')
plt.xscale('log')
plt.yscale('log')
plt.xlim(5e-3, 1.25)
plt.ylim(1e-2, 400)
plt.title('density')


P = (gamma - 1.0) * rho * u
Px = (gamma - 1.0) * rhox * ux

plt.subplot(2, 3, 3)
plt.scatter(rr, P, s = 0.05, color = 'orange')
plt.scatter(rx, Px, s = 20, color = 'b')
plt.plot(radius_ppm, press, linestyle = '--', color = 'k')
plt.xscale('log')
plt.yscale('log')
plt.xlim(5e-3, 1.25)
plt.ylim(1e-2, 400)
plt.title('pressure')


plt.subplot(2, 3, 4)
plt.scatter(rr, u, s = 0.05, color = 'orange')
plt.scatter(rx, ux, s = 20.0, color = 'b')
plt.plot(radius_ppm, u_ppm, linestyle = '--', color = 'k')
plt.xscale('log')
plt.yscale('log')
plt.xlim(5e-3, 1.25)
plt.ylim(1e-2, 3)
plt.title('internal energy')


ent_s = (gamma - 1.0) * rho**(1.0 - gamma) * u
ent_sx = (gamma - 1.0) * rhox**(1.0 - gamma) * ux

plt.subplot(2, 3, 5)
plt.scatter(rr, ent_s, s = 0.05, color = 'orange')
plt.scatter(rx, ent_sx, s = 20.0, color = 'b')
plt.plot(radius_ppm, entropy, linestyle = '--', color = 'k')
plt.xscale('log')
plt.xlim(5e-3, 1.25)
plt.ylim(0, 0.25)
plt.title('entropy')


plt.savefig('v_u_rho_vs_r.png')

plt.show()




