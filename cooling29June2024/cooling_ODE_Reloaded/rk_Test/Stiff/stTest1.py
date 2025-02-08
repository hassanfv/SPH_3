import numpy as np
import matplotlib.pyplot as plt

# Define the system of differential equations
def system(x, U):
    u, v = U
    du_dx = 998*u + 1998*v
    dv_dx = -999*u - 1999*v
    return np.array([du_dx, dv_dx])

# Implementing the 4th-order Runge-Kutta method
def runge_kutta(f, U0, x_range, h):
    x_values = np.arange(x_range[0], x_range[1], h)
    U_values = np.zeros((len(x_values), len(U0)))
    
    U = np.array(U0)
    for i, x in enumerate(x_values):
        U_values[i] = U
        k1 = f(x, U)
        k2 = f(x + h/2, U + h*k1/2)
        k3 = f(x + h/2, U + h*k2/2)
        k4 = f(x + h, U + h*k3)
        U = U + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
    
    return x_values, U_values

# Exact solution
def exact_solution(x):
    u_exact = 2*np.exp(-x) - np.exp(-1000*x)
    v_exact = -np.exp(-x) + np.exp(-1000*x)
    return u_exact, v_exact

# Initial conditions
U0 = [1, 0]
x_range = (0, 0.05)  # Keeping a short range due to stiffness
h = 0.0001  # Small step size to capture stiffness

# Solve using Runge-Kutta
x_values, U_values = runge_kutta(system, U0, x_range, h)
u_rk, v_rk = U_values[:, 0], U_values[:, 1]

# Compute exact solutions
u_exact_values, v_exact_values = exact_solution(x_values)

# Plot results
plt.figure(figsize=(10, 5))
plt.plot(x_values, u_rk, label='Runge-Kutta u(x)', linestyle='dashed')
plt.plot(x_values, v_rk, label='Runge-Kutta v(x)', linestyle='dashed')
plt.plot(x_values, u_exact_values, label='Exact u(x)', linewidth=2)
plt.plot(x_values, v_exact_values, label='Exact v(x)', linewidth=2)
plt.xlabel('x')
plt.ylabel('Solution values')
plt.legend()
plt.title('Comparison of Runge-Kutta and Exact Solution')
plt.show()

