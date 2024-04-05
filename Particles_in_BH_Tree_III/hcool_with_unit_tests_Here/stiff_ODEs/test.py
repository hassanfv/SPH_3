import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Define the rate constants
k1 = 0.04
k2 = 3e7
k3 = 1e4

# Define the system of ODEs
def robertson(t, y):
    A, B, C = y
    dAdt = -k1 * A + k3 * B * C
    dBdt = k1 * A - k2 * B**2 - k3 * B * C
    dCdt = k2 * B**2
    return [dAdt, dBdt, dCdt]

# Initial conditions
y0 = [1.0, 0.0, 0.0]

# Time span for the solution
t_span = (0, 40)

# Solve the system of ODEs
solution = solve_ivp(robertson, t_span, y0, method='LSODA', dense_output=True)

# Plot the solution
t = np.linspace(t_span[0], t_span[1], 400)
y = solution.sol(t)

plt.plot(t, y.T)
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.title('Robertson Problem Solution')
plt.legend(['A', 'B', 'C'])
plt.grid()
plt.show()

