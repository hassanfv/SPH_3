import numpy as np
import matplotlib.pyplot as plt

# Parameters
h = 0.2  # Time step
t_max = 4  # Maximum time
N = int(t_max / h)  # Number of steps
tol = 1e-10  # Tolerance for Newton-Raphson
max_iter = 1000  # Maximum iterations for Newton-Raphson

# Initial condition
t = np.linspace(0, t_max, N + 1)
x = np.zeros(N + 1)
x[0] = 1  # Initial value

# Function definitions
def f(t, x):
    return (1 - 2 * t) * x

def df_dx(t, x):
    return 1 - 2 * t

# Newton-Raphson method to solve the implicit equation
def newton_raphson(t_n1, x_n, h):
    x_guess = x_n  # Initial guess
    for i in range(max_iter):
        F = x_guess - x_n - h * f(t_n1, x_guess)
        J = 1 - h * df_dx(t_n1, x_guess)  # Jacobian (df/dx)
        delta = F / J
        x_guess -= delta
        if abs(delta) < tol:
            break
    else:
        raise RuntimeError("Newton-Raphson did not converge")
    return x_guess

# Backward Euler method
for n in range(N):
    t_n1 = t[n + 1]
    x_n = x[n]
    x[n + 1] = newton_raphson(t_n1, x_n, h)

# Exact solution
x_exact = np.exp(1./4. - (t - 1./2.)**2)

# Plotting the results
plt.figure(figsize=(10, 6))
plt.plot(t, x, 'bo-', label='Backward Euler Solution')
plt.plot(t, x_exact, 'r--', label='Exact Solution')
plt.xlabel('t')
plt.ylabel('x(t)')
plt.title('Comparison of Backward Euler Solution and Exact Solution')
plt.legend()
plt.grid(True)
plt.show()

