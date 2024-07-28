
import numpy as np
from scipy.linalg import solve

# Define the function f(u, t) representing the system of ODEs
def f(u, t):
    return np.array([
        -u[0] + 2*u[1],
        -3*u[1] + u[2],
        4*u[2] - 2*u[3],
        -u[3] + u[4],
        -2*u[4] + u[0]
    ])

# Define the function g(u, u_n, dt, t) based on implicit Euler
def g(u, u_n, dt, t):
    return u - u_n - dt * f(u, t)

# Compute the Jacobian J of g with respect to u
def jacobian(u, dt):
    n = len(u)
    J = np.zeros((n, n))
    eps = 1e-8
    for i in range(n):
        for j in range(n):
            u1 = u.copy()
            u1[j] += eps
            J[i, j] = (g(u1, u, dt, t)[i] - g(u, u, dt, t)[i]) / eps
    return J

# Newton's method to solve the system
def newton_method(u_n, dt, t, tol=1e-8, max_iter=100):
    u_k = u_n
    for _ in range(max_iter):
        J = jacobian(u_k, dt)
        g_val = g(u_k, u_n, dt, t)
        a = solve(J, g_val)  # Solve Ja = g(x_k) for a
        u_k = u_k - a  # Update u_k
        if np.linalg.norm(a) < tol:
            break
    return u_k

# Initial conditions and parameters
u_n = np.array([1.0, 0.5, -0.5, 1.0, -1.0])
dt = 0.1
t = 0.0

# Solve for u_n+1 using Newton's method
u_n1 = newton_method(u_n, dt, t)
print("Solution u_n+1:", u_n1)

print()
print()
print('------ Test -------')
JJ = np.array([[1, 2, 1], [0, 3, 2], [1, 4, 3]])
gg = np.array([8, 12, 18])
print(solve(JJ, gg))

print()
print()
print('------ Test2 -------')
JJ = np.array([[1, 2], [2, 3]])
gg = np.array([4, 7])
print(solve(JJ, gg))



