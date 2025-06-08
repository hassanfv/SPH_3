
# Only Forward Euler

import numpy as np
import matplotlib.pyplot as plt


def func(t):
  return -2.0 * 3.0 * np.exp(-2. * t)

xn = 3.0
tn = 0.0

res = [[tn, xn]]

h = 0.10

while tn < 2:

  xn = xn + h * func(tn)
  tn += h

  res.append([tn , xn])

res = np.array(res)

x = np.linspace(0, 2, 100)
y = 3.0 * np.exp(-2. * x)
plt.plot(x, y, color = 'b', label = "Analytical")

#yy = 3.0 - 6 * x
#plt.plot(x, yy, color = 'r')

plt.plot(res[:, 0], res[:, 1], label = "Numerical", color = 'limegreen', marker = 'o')

#plt.ylim(-1, 4)

plt.title(f"h = {h}")

plt.legend()

plt.show()

