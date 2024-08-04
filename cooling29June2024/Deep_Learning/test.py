
import numpy as np
import pandas as pd
import pickle


with open('Test.pkl', 'rb') as f:\
  data = pickle.load(f)


inRes = data['inRes']
outRes = data['outRes']
y_pred = data['y_pred']

#print(inRes.shape, outRes.shape, y_pred.shape)

i = 0

T_before = inRes[i, 1]

T_after_real = outRes[i, 1]

T_after_pred = y_pred[i, 1]

print(f'T_before = {T_before},  T_after_real = {T_after_real},  T_after_pred = {T_after_pred}')

