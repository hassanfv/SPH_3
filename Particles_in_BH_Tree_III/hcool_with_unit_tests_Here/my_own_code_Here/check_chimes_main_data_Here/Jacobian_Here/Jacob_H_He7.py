#import h5py
import numpy as np
import matplotlib.pyplot as plt
#from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp
import pickle
from data import *


#----- dT_dt
def dT_dt(nH0, nHp, nHe0, nHep, nHepp, T):
    Tx = np.log10(T)
    ne = nHp + nHep + 2.0 * nHepp
    ntot = nH0 + nHp + nHe0 + nHep + nHepp + ne
    Lamb = (
        10**g1(Tx) * ne * nH0  # H0
        + 10**g2(Tx) * ne * nHp # Hp
        + 10**g3(Tx) * nHe0 * ne # He0 
        + 10**g4(Tx) * nHep * ne # Hep 
        + 10**g5(Tx) * nHepp * ne# Hepp
    )
    dT_dt = -1.0 * (gamma - 1.0) / kB / ntot * Lamb
    return dT_dt
  
#----- gradient function
def gradient(f, y, idx, delta):
    y_up = y.copy()
    y_up[idx] += delta
    return (f(*y_up) - f(*y)) / delta

#----- Jacobian function
def jacobian(t, y):
    nH0, nHp, nHe0, nHep, nHepp, T = y
    delta = 1e-6

    # List of functions corresponding to each differential equation
    funcs = [dnH0_dt, dnHp_dt, dnHe0_dt, dnHep_dt, dnHepp_dt, dT_dt]

    # Initialize Jacobian matrix
    jacobian_matrix = np.zeros((6, 6))

    for i in range(6):
        for j in range(6):
            jacobian_matrix[i, j] = gradient(funcs[i], y, j, delta)
    
    return jacobian_matrix

#----- dnH0_dt
def dnH0_dt(nH0, nHp, nHe0, nHep, nHepp, T):
    Tx = np.log10(T)
    ne = nHp + nHep + 2.0 * nHepp
    return 10**k2(Tx) * nHp * ne - 10**k1(Tx) * nH0 * ne

#----- dnHp_dt
def dnHp_dt(nH0, nHp, nHe0, nHep, nHepp, T):
    Tx = np.log10(T)
    ne = nHp + nHep + 2.0 * nHepp
    return 10**k1(Tx) * nH0 * ne - 10**k2(Tx) * nHp * ne

#----- dnHe0_dt
def dnHe0_dt(nH0, nHp, nHe0, nHep, nHepp, T):
    Tx = np.log10(T)
    ne = nHp + nHep + 2.0 * nHepp
    return 10**k5(Tx) * nHep * ne - 10**k3(Tx) * nHe0 * ne

#----- dnHep_dt
def dnHep_dt(nH0, nHp, nHe0, nHep, nHepp, T):
    Tx = np.log10(T)
    ne = nHp + nHep + 2.0 * nHepp
    return 10**k6(Tx) * nHepp * ne + 10**k3(Tx) * nHe0 * ne - 10**k4(Tx) * nHep * ne - 10**k5(Tx) * nHep * ne

#----- dnHepp_dt
def dnHepp_dt(nH0, nHp, nHe0, nHep, nHepp, T):
    Tx = np.log10(T)
    ne = nHp + nHep + 2.0 * nHepp
    return 10**k4(Tx) * nHep * ne - 10**k6(Tx) * nHepp * ne

#----- Lambda
def Lambda(T, nH0, nHp, nHe0, nHep, nHepp):
    Tx = np.log10(T)
    ne = nHp + nHep + 2.0 * nHepp
    Lamb = (
        10**g1(Tx) * ne * nH0  # H0
        + 10**g2(Tx) * ne * nHp # Hp
        + 10**g3(Tx) * nHe0 * ne # He0 
        + 10**g4(Tx) * nHep * ne # Hep 
        + 10**g5(Tx) * nHepp * ne# Hepp
    )
    return Lamb

#----- func
def func(t, y):
    nH0, nHp, nHe0, nHep, nHepp, T = y
    Tx = np.log10(T)
    ne = nHp + nHep + 2.0 * nHepp
    ntot = nH0 + nHp + nHe0 + nHep + nHepp + ne
    dnH0_dt = 10**k2(Tx) * nHp * ne - 10**k1(Tx) * nH0 * ne
    dnHp_dt = 10**k1(Tx) * nH0 * ne - 10**k2(Tx) * nHp * ne
    dnHe0_dt = 10**k5(Tx) * nHep * ne - 10**k3(Tx) * nHe0 * ne
    dnHep_dt = 10**k6(Tx) * nHepp * ne + 10**k3(Tx) * nHe0 * ne - 10**k4(Tx) * nHep * ne - 10**k5(Tx) * nHep * ne
    dnHepp_dt = 10**k4(Tx) * nHep * ne - 10**k6(Tx) * nHepp * ne
    Lamb = Lambda(T, nH0, nHp, nHe0, nHep, nHepp)
    dT_dt = -1.0 * (gamma - 1.0) / kB / ntot * Lamb
    return [dnH0_dt, dnHp_dt, dnHe0_dt, dnHep_dt, dnHepp_dt, dT_dt]

nH = 1000.0
He_solar = 10**(-1.07)
nHe = He_solar * nH
print('nHe (cm^-3) = ', nHe)

nH0_i = 1.0
nHp_i = 999.0
nHe0_i = 0.01 * nHe
nHep_i = 0.1 * nHe
nHepp_i = nHe - nHe0_i - nHep_i

y0 = [nH0_i, nHp_i, nHe0_i, nHep_i, nHepp_i, 1e6]
t_span = (1*3.16e7, 5000*3.16e7)

# Solve with Jacobian
solution_with_jac = solve_ivp(func, t_span, y0, method='LSODA', jac=jacobian, dense_output=True)

# Solve without Jacobian
solution_without_jac = solve_ivp(func, t_span, y0, method='LSODA', dense_output=True)

# Define the time points where you want the solution
t = np.linspace(t_span[0], t_span[1], 10000)

# Evaluate the solutions at the defined time points
y_with_jac = solution_with_jac.sol(t)
y_without_jac = solution_without_jac.sol(t)

# Calculate the differences between the solutions
differences = np.abs(y_with_jac - y_without_jac)

# Performance comparison
print("Performance with Jacobian:")
print("Number of function evaluations:", solution_with_jac.nfev)
print("Number of Jacobian evaluations:", solution_with_jac.njev)
print("Number of LU decompositions:", solution_with_jac.nlu)

print("\nPerformance without Jacobian:")
print("Number of function evaluations:", solution_without_jac.nfev)
print("Number of Jacobian evaluations:", solution_without_jac.njev)
print("Number of LU decompositions:", solution_without_jac.nlu)


res = y_with_jac.T

print(res.shape)

t_yrs = t / 3.16e7

T = res[:, -1]

nH0 = res[:, 0]
nHp = res[:, 1]
nHe0 = res[:, 2]
nHep = res[:, 3]
nHepp= res[:, 4]

#------ Result from "test_primordial_hdf5_v2.py" code -----
with open('chimesRes_H.pkl', 'rb') as f:
  df = pickle.load(f)
# dictx = {'t_Arr_in_yrs': t_Arr_in_yrs, 'TEvol': TEvol, 'nHe0': nHe0, 'nHep': nHep, 'nHepp': nHepp}
t_Arr_in_yrsx = df['t_Arr_in_yrs']
TEvolx = df['TEvol']
nH0x = df['nH0']
nHpx = df['nHp']
nHe0x = df['nHe0']
nHepx = df['nHep']
nHeppx = df['nHepp']
nHeTotx = nHe0x + nHepx + nHeppx
#----------------------------------------------------------


plt.figure(figsize = (16, 8))

plt.subplot(2, 3, 1)
plt.scatter(t_yrs, np.log10(T), s = 5, color = 'k', label = 'my own code')
plt.scatter(t_Arr_in_yrsx, np.log10(TEvolx), s = 2, color = 'orange', label = 'chimes result', linestyle = '--')
plt.xlim(0, 3000)
plt.ylim(3, 8)
plt.legend()

plt.subplot(2, 3, 2)

plt.plot(t_yrs, nHe0/nHe, color = 'r', label = 'nHe0')
plt.plot(t_yrs, nHep/nHe, color = 'g', label = 'nHep')
plt.plot(t_yrs, nHepp/nHe, color = 'b', label = 'nHepp')

plt.plot(t_Arr_in_yrsx, nHe0x/nHeTotx, color = 'r', label = 'nHe0 - chimes', linestyle = ':')
plt.plot(t_Arr_in_yrsx, nHepx/nHeTotx, color = 'g', label = 'nHep - chimes', linestyle = ':')
plt.plot(t_Arr_in_yrsx, nHeppx/nHeTotx, color = 'b', label = 'nHepp - chimes', linestyle = ':')

plt.xlim(0, 5000)
plt.ylim(1e-7, 2.)

plt.yscale('log')
plt.title('solve_ivp')
plt.legend()


plt.subplot(2, 3, 3)
nHeTot = nHe0 + nHep + nHepp
plt.plot(T, nH0/nH, label = 'nH0', color = 'r')
plt.plot(T, nHp/nH, label = 'nHp', color = 'g')
plt.plot(T, nHe0/nHeTot, label = 'nHe0', color = 'b')
plt.plot(T, nHep/nHeTot, label = 'nHep', color = 'orange')
plt.plot(T, nHepp/nHeTot,label = 'nHepp', color = 'purple')

plt.plot(TEvolx, nH0x/nH, label = 'nH0 - chimes', color = 'r', linestyle = ':')
plt.plot(TEvolx, nHpx/nH, label = 'nHp - chimes', color = 'g', linestyle = ':')
plt.plot(TEvolx, nHe0x/nHeTotx, label = 'nHe0 - chimes', color = 'b', linestyle = ':')
plt.plot(TEvolx, nHepx/nHeTotx, label = 'nHep - chimes', color = 'orange', linestyle = ':')
plt.plot(TEvolx, nHeppx/nHeTotx,label = 'nHepp - chimes', color = 'purple', linestyle = ':')

plt.yscale('log')
plt.xscale('log')
plt.ylim(2e-3, 1.2)
plt.xlim(1e4, 5e6)
plt.legend()


plt.tight_layout()

plt.savefig('HeH_jacob.png')

plt.show()





