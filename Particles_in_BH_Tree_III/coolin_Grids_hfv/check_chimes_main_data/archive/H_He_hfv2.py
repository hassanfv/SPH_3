import h5py
import numpy as np
import matplotlib.pyplot as plt


with h5py.File('chimes_main_data.hdf5', 'r') as file:
  rates = file['T_dependent/rates'][:]
  ratesAB = file['recombination_AB/rates_caseA'][:]
  Temp = file['TableBins/Temperatures'][:]
  cooling_rates = file['cooling/rates'][:]

# NOTE: I used "T_dependent_reactants.py" code to find indices 111, 0, 108, etc!!!

#---- Reaction rates -------
k1x = rates[111, :] # Reaction: (H0 + e ---> Hp + 2e) ::: H0 Collisional ionization
k2x = ratesAB[0, :] # Reaction: (Hp + e ---> H0 + γ)  ::: photo-recombination of Hp and e.
k3x = rates[108, :] # Reaction: (He0 + e ---> Hep + 2e) ::: He0 Collisional ionization 
k4x = rates[0, :]   # Reaction: (Hep + e ---> Hepp + 2e) ::: Hep Collisional ionization
k5x = ratesAB[1, :] # Reaction: (Hep + e ---> He0 + γ)  ::: photo-recombination of Hep and e. !!! PROBABLY di-electric is also included ?????
k6x = rates[222, :] # Reaction: (Hepp + e ---> Hep + γ)  ::: photo-recombination of Hepp and e.


#---- Cooling rates ------
g1x = cooling_rates[0, :] # cooling via H0
g2x = cooling_rates[1, :] # cooling via Hp

g3x = cooling_rates[2, :] # cooling via He0
g4x = cooling_rates[3, :] # cooling via Hep
g5x = cooling_rates[4, :] # cooling via Hepp

plt.figure(figsize = (14, 6))

plt.subplot(1, 2, 1)
plt.scatter(Temp, k1x, s = 2, label = 'H0 + e ---> Hp + 2e')
plt.scatter(Temp, k2x, s = 2, label = 'Hp + e ---> H0 + γ')
plt.scatter(Temp, k3x, s = 2, label = 'He0 + e ---> Hep + 2e')
plt.scatter(Temp, k4x, s = 2, label = 'Hep + e ---> Hepp + 2e')
plt.scatter(Temp, k5x, s = 2, label = 'Hep + e ---> He0 + γ')
plt.scatter(Temp, k6x, s = 2, label = 'Hepp + e ---> Hep + γ')
plt.ylim(-26, -6)
plt.legend()


plt.subplot(1, 2, 2)
plt.scatter(Temp, g1x, s = 2, label = 'cooling via H0')
plt.scatter(Temp, g2x, s = 2, label = 'cooling via Hp')

plt.scatter(Temp, g3x, s = 2, label = 'cooling via He0')
plt.scatter(Temp, g4x, s = 2, label = 'cooling via Hep')
plt.scatter(Temp, g5x, s = 2, label = 'cooling via Hepp')
plt.legend()

plt.show()


