
import numpy as np
import matplotlib.pyplot as plt
import pickle
import pandas as pd
from scipy.interpolate import CubicSpline
from numpy.polynomial.polynomial import Polynomial
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import interp1d


with open('test.pkl', 'rb') as f:
  data = pickle.load(f)

t_yrs = data['t_Arr_in_yrs']
TEvol = data['TEvol']


df = pd.read_csv('data.csv')
t = df['t']
T = df['T']

# Linear interpolation
linear_interp = interp1d(t, T, kind='linear')

# Generate fine-grained t values for smooth plotting
t_fine = np.linspace(t.min(), t.max(), 500)
T_fine = linear_interp(t_fine)


#plt.scatter(t_yrs, TEvol, s = 1, color = 'k')
plt.scatter(t, T, color = 'red', s = 10)
plt.scatter(t_fine, T_fine, color = 'k', s = 5)
plt.show()


