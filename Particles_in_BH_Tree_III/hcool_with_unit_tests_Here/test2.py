
import numpy as np
from photolibs3 import *
import matplotlib.pyplot as plt
import pandas as pd



T_n = 1e5

nH = 100.0

Tgrid = np.logspace(4.0, 9.0, 100)

res = []

for T_n in Tgrid:
  Gam, Lamb_at_T_n = coolingHeatingRates(T_n, nH)

  res.append(Lamb_at_T_n)

res = np.array(res)

#-------------- From Gnat & Sternberg - 2007 ---------------------
df = pd.read_csv('table4a.dat', delim_whitespace=True, header=None)
T = df.iloc[:, 0].values
Lam1 = df.iloc[:, 1].values
Lam3 = df.iloc[:, 3].values
#------------------------------------

plt.scatter(np.log10(Tgrid), np.log10(res), s = 10)
plt.plot(np.log10(T), np.log10(Lam1), color = 'k')
plt.plot(np.log10(T), np.log10(Lam3), color = 'orange')

plt.show()





