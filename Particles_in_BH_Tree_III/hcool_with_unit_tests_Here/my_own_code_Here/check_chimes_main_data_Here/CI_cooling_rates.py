import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import interp1d


#-----------------------------------------------------------
#------------> Cooling rates from CHIMES table <------------
#-----------------------------------------------------------
with h5py.File('chimes_main_data.hdf5', 'r') as file:

  coolants_4d = file['cooling/coolants_4d'][:] # is the same for low and hiT.... Only Temp will change!! And in hiT it is only a function of T!!!
  
  Temp_4d = file['TableBins/cool_4d_Temperatures'][:]
  HIDensity_4d = file['TableBins/cool_4d_HIDensities'][:]
  elecDensity_4d = file['TableBins/cool_4d_ElectronDensities'][:]
  HIIDensity_4d = file['TableBins/cool_4d_HIIDensities'][:]
  
  rates_4d = file['cooling/rates_4d'][:] # NOTE it is rates_4d for low Temp
  
  #-------- hiT_4d ---------------
  Temp_hiT_4d = file['TableBins/cool_hiT_4d_Temperatures'][:]
  rates_hiT_4d = file['cooling/rates_hiT_4d'][:] # NOTE it is rates_4d for high Temp
  
  
  print(elecDensity_4d)
  print()
  print(HIDensity_4d)
  print()
  print(HIIDensity_4d)
  print()
  
  print(Temp_4d.shape, rates_4d.shape, Temp_hiT_4d.shape, rates_hiT_4d.shape)


print()
print('coolants_4d.shape: ', coolants_4d.shape)
print()
print('rates_4d.shape: ', rates_4d.shape)
print()
print('coolants_4d = ', coolants_4d)
print()
print('rates_hiT_4d.shape: ', rates_4d.shape)
print()

nt = np.where(coolants_4d == 7)[0]   # Change the index here to see the reactions and products of that element
print('******** coolants ****************')
print('nt = ', nt)
print(coolants_4d[nt])
print('************************')
print()


print('Temp_4d = ', Temp_4d)
print()

C0_rate_4d = rates_4d[0, :]
interp_4d = RegularGridInterpolator((Temp_4d, HIDensity_4d, elecDensity_4d, HIIDensity_4d), C0_rate_4d)

test_lowT = np.array([3.0, -8.0, -8.0, -8.0])
print('lowT test = ', interp_4d(test_lowT))


C0_rate_hiT_4d = rates_hiT_4d[0, :]
interp_hiT_4d = interp1d(Temp_hiT_4d, C0_rate_hiT_4d, kind='linear', fill_value="extrapolate")

test_hiT_T = 5.15
print('hiT test = ', interp_hiT_4d(test_hiT_T))



plt.scatter(Temp_4d, rates_4d[0, :, 0, 0, 0], s = 2, label = 'C0 - low_T_4d')
plt.scatter(Temp_hiT_4d, rates_hiT_4d[0, :], s = 2, label = 'C0 - hi_T_4d')

plt.legend()
plt.show()




