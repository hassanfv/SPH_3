
# Forward Euler and RK2

import numpy as np
import matplotlib.pyplot as plt


def func(t):
  return -2.0 * 3.0 * np.exp(-2. * t)

xn = xnRK2 = 3.0
tn = 0.0

res = [[tn, xn, xnRK2]]

h = 0.10

while tn < 2:

  xn = xn + h * func(tn)
  xnRK2 = xnRK2 + h * func(tn + h/2.)
  tn += h

  res.append([tn , xn, xnRK2])

res = np.array(res)

x = np.linspace(0, 2, 100)
y = 3.0 * np.exp(-2. * x)
plt.plot(x, y, color = 'b', label = "Analytical", linewidth = 3)

plt.plot(res[:, 0], res[:, 1], label = "F-Euler", color = 'limegreen', marker = 'o')
plt.plot(res[:, 0], res[:, 2], label = "RK2", color = 'crimson', marker = 'o')

plt.title(f"h = {h}")

plt.legend()

plt.show()

