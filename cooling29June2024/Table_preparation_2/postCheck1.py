
import numpy as np
#import pandas as pd
import pickle
import matplotlib.pyplot as plt


with open('TEvol.pkl', 'rb') as f:
  data = pickle.load(f)

#['TEvol', 'nH', 'Temp', 'rkpc', 'Lsh', 't']

TEvol = data['TEvol']
nH = data['nH']
Temp = data['Temp']
rkpc = data['rkpc']
Lsh = data['Lsh'] # ---- This is in log10 of parsec ! So 0 means 1.0 pc!
t = data['t']

print(TEvol.shape)

print()
print(f'nH = {nH}\n')
print(f'Temp = {Temp}\n')
print(f'rkpc = {rkpc}\n')
print(f'Lsh = {Lsh}\n')


