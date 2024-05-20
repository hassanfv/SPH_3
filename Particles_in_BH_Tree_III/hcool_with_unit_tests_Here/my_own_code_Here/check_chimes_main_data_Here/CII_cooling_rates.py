import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import interp1d


#-----------------------------------------------------------
#------------> Cooling rates from CHIMES table <------------
#-----------------------------------------------------------
with h5py.File('chimes_main_data.hdf5', 'r') as file:

  coolants_2d = file['cooling/coolants_2d'][:] # is the same for low and hiT.... Only Temp will change!! And in hiT it is only a function of T!!!
  
  Temp_2d = file['TableBins/cool_2d_Temperatures'][:]
  elecDensity_2d = file['TableBins/cool_2d_ElectronDensities'][:]
  
  rates_2d = file['cooling/rates_2d'][:] # NOTE it is rates_2d for low Temp
  
  #-------- hiT_2d ---------------
  Temp_hiT_2d = file['TableBins/cool_hiT_2d_Temperatures'][:]
  rates_hiT_2d = file['cooling/rates_hiT_2d'][:] # NOTE it is rates_2d for high Temp
  
  
  print(elecDensity_2d)
  print()
  
  print(Temp_2d.shape, rates_2d.shape, Temp_hiT_2d.shape, rates_hiT_2d.shape)


print()
print('coolants_2d.shape: ', coolants_2d.shape)
print()
print('rates_2d.shape: ', rates_2d.shape)
print()
print('coolants_2d = ', coolants_2d)
print()
print('rates_hiT_2d.shape: ', rates_2d.shape)
print()


nt = np.where(coolants_2d == 8)[0]   # 8 is C+, Change the index here to see the reactions and products of that element
print('******** coolants ****************')
print('nt = ', nt)
print(coolants_2d[nt])
print('************************')
print()


print('Temp_2d = ', Temp_2d)
print()

Cp_rate_2d = rates_2d[0, :]
interp_2d = RegularGridInterpolator((Temp_2d, elecDensity_2d), Cp_rate_2d)

test_lowT = np.array([3.0, -8.0])
print('lowT test = ', interp_2d(test_lowT))


Cp_rate_hiT_2d = rates_hiT_2d[0, :]
interp_hiT_2d = interp1d(Temp_hiT_2d, Cp_rate_hiT_2d, kind='linear', fill_value="extrapolate")

test_hiT_T = 5.15
print('hiT test = ', interp_hiT_2d(test_hiT_T))



plt.scatter(Temp_2d, rates_2d[0, :, 0], s = 2, label = 'Cp - low_T_2d')
plt.scatter(Temp_hiT_2d, rates_hiT_2d[0, :], s = 2, label = 'Cp - hi_T_2d')

plt.legend()
plt.show()




