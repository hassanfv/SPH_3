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
t_span = (1*3.16e7, 20000*3.16e7)

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

# Plot the results for comparison
plt.figure(figsize=(18, 6))

# Plot y0, y1, y2 with Jacobian
plt.subplot(1, 3, 1)
plt.plot(t, y_with_jac.T)
plt.xscale('log')
plt.title('Solution with Jacobian')
plt.xlabel('Time')
plt.ylabel('Concentrations')
plt.legend(['nH0', 'nHp', 'nHe0', 'nHep', 'nHepp', 'T'])

# Plot y0, y1, y2 without Jacobian
plt.subplot(1, 3, 2)
plt.plot(t, y_without_jac.T)
plt.xscale('log')
plt.title('Solution without Jacobian')
plt.xlabel('Time')
plt.ylabel('Concentrations')
plt.legend(['nH0', 'nHp', 'nHe0', 'nHep', 'nHepp', 'T'])

# Plot differences
plt.subplot(1, 3, 3)
plt.plot(t, differences.T)
plt.xscale('log')
plt.yscale('log')
plt.title('Differences between solutions')
plt.xlabel('Time')
plt.ylabel('Absolute Difference')
plt.legend(['nH0', 'nHp', 'nHe0', 'nHep', 'nHepp', 'T'])

plt.tight_layout()
plt.show()


