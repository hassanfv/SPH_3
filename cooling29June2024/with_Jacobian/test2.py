import numpy as np

# Define constants
N = 5  # Number of ODEs
M = 4  # Number of independent variables
H = 1e-5  # Step size for finite difference

def ode_system(x):
    # Define your system of ODEs here
    # Example:
    f = np.zeros(N)
    f[0] = x[0] * x[1] + x[2] * x[3]
    f[1] = x[1] * x[1] - x[0]
    f[2] = x[2] + x[0] * x[3]
    f[3] = x[3] - x[1] * x[2]
    f[4] = x[0] + x[1] + x[2] + x[3]
    return f

def compute_jacobian(x):
    J = np.zeros((N, M))
    for col in range(M):
        x_plus = np.copy(x)
        x_minus = np.copy(x)
        
        x_plus[col] += H
        x_minus[col] -= H

        f_plus = ode_system(x_plus)
        f_minus = ode_system(x_minus)

        for row in range(N):
            J[row, col] = (f_plus[row] - f_minus[row]) / (2.0 * H)
    return J


x = np.array([1.0, 2.0, 3.0, 4.0])  # Example initial point
J = compute_jacobian(x)

# Print the Jacobian matrix
print("Jacobian matrix:")
for i in range(N):
    for j in range(M):
        print(f"{J[i, j]: .5f}", end=" ")
    print()




