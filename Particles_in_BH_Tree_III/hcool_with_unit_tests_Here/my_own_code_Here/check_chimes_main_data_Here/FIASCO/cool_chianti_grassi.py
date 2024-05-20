import os
os.environ['XUVTOP'] = '/home/pc/CHIANTI/'

import ChiantiPy
import ChiantiPy.core as ch
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

# temperature range
trange = np.logspace(4, 8., 15)

# loop on ions, 1=neutral
for ll in tqdm(range(1, 7)):

    ems = np.zeros_like(trange)  # store emission here
    pdens = np.ones_like(trange)  # proton density

    # loop on temperature
    for i, t in enumerate(trange):
        # create atom and compute emission, CHIANTI carbon netural is c_1
        h1 = ch.ion("c_%d" % ll, temperature=t, eDensity=1e0, abundance=1e0, pDensity=pdens)
        h1.emiss()
        em = h1.Emiss["emiss"] * 4e0 * np.pi  # CHIANTI is in erg/s/cm3/sr, this is erg/s/cm3
        ems[i] = np.clip(np.sum(em, axis=0), 1e-40, 1e99)  # sums over all the emissions

    # plot using neutral=0 labelling as in Gnat paper
    plt.loglog(trange, ems, label=ll-1)

plt.ylim(1e-24, 1e-17)

plt.legend()
plt.show()
