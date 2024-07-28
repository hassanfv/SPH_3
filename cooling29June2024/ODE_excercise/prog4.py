
# Implicit

import numpy as np
import matplotlib.pyplot as plt


def func(t_n, x_n):
  return (1. - 2. * t_n) * x_n 

def g(x):
  return x - xn - h * func(tn, xn)

def deriv_g(x):
  return 1. - h * (1. - 2 * (tn))


def find_root(x):

  max_iterations = 1000
  for i in range(max_iterations):
    fx = g(x)
    dfx = deriv_g(x)
    
    # Update x using the Newton-Raphson formula
    x_new = x - fx / dfx
    
    x = x_new

  return x

tn = 0.
xn = 1.

h = 0.1

res = []
res.append([tn, xn])

tn +=h

while tn <= 4.:
  
  xnp1 = find_root(xn)
  
  res.append([tn, xnp1])
  
  xn = xnp1
  tn += h

res = np.array(res)

tt = res[:, 0]
xx = res[:, 1]

t = np.arange(0, 4, 0.001)
x = np.exp(1./4. - (t - 1./2.)**2)

plt.plot(t, x, color = 'k')
plt.scatter(tt, xx, color = 'orange')


plt.show()
  
  

