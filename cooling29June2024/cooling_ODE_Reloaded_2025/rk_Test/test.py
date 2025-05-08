import numpy as np
import matplotlib.pyplot as plt

# Define the differential equation dy/dt = f(t, y)
def f(t, y):
    return -y

# Implement the Runge-Kutta 4th order method
def runge_kutta(f, t0, y0, t_end, h):
    t_values = np.arange(t0, t_end + h, h)
    y_values = np.zeros(len(t_values))
    y_values[0] = y0
    
    for i in range(1, len(t_values)):
        t = t_values[i - 1]
        y = y_values[i - 1]
        
        k1 = h * f(t, y)
        k2 = h * f(t + h / 2, y + k1 / 2)
        k3 = h * f(t + h / 2, y + k2 / 2)
        k4 = h * f(t + h, y + k3)
        
        y_values[i] = y + (k1 + 2 * k2 + 2 * k3 + k4) / 6
    
    return t_values, y_values

# Implement the Forward Euler method
def forward_euler(f, t0, y0, t_end, h):
    t_values = np.arange(t0, t_end + h, h)
    y_values = np.zeros(len(t_values))
    y_values[0] = y0
    
    for i in range(1, len(t_values)):
        t = t_values[i - 1]
        y = y_values[i - 1]
        
        y_values[i] = y + h * f(t, y)
    
    return t_values, y_values

# Parameters
t0 = 0
y0 = 1
t_end = 5
h = 0.1

# Solve using Runge-Kutta and Forward Euler
t_rk, y_rk = runge_kutta(f, t0, y0, t_end, h)
t_fe, y_fe = forward_euler(f, t0, y0, t_end, 0.1*h)

# Analytical solution
t_exact = np.linspace(t0, t_end, 1000)
y_exact = np.exp(-t_exact)

# Plot results
plt.figure(figsize=(10, 6))
plt.plot(t_exact, y_exact, label="Exact Solution", color="black", linestyle="--")
plt.plot(t_rk, y_rk, label="Runge-Kutta 4th Order", marker="o")
plt.plot(t_fe, y_fe, label="Forward Euler", marker="x")
plt.xlabel("Time (t)")
plt.ylabel("y(t)")
plt.title("Comparison of Runge-Kutta and Forward Euler Methods")
plt.legend()
plt.grid()
plt.show()

