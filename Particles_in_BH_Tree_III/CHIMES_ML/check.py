
import numpy as np
import pickle
import pandas as pd

#--- InputList ------
with open('inputListsX.pkl', 'rb') as f:
  inputData = pickle.load(f)

print('inputData.shape (T, nH, rkpc, NH) = ', inputData.shape)
print()


df = pd.DataFrame(inputData)
df.to_csv('data.csv', index = False)



with open('CHIMES_EC2_hfv.pkl', 'rb') as f:
  dictx = pickle.load(f)

# ['Tres', 'uRes', 'muRes', 'AbRes']

T = np.array(dictx['Tres'])
u = np.array(dictx['uRes'])
mu = np.array(dictx['muRes'])
Ab = np.array(dictx['AbRes'])

print('T.shape = ', T.shape)
print('u.shape = ', u.shape)
print('mu.shape = ', mu.shape)
print('Ab.shape = ', Ab.shape)
print()


j = 95
print('input: ', inputData[j, :])
print()
print(T[j, :])


