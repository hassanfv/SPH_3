
import numpy as np
import tensorflow as tf
import pickle
import time


with open ("initialData.pkl", "rb") as f:
  data = pickle.load(f)

nn = np.where(data[:, 0] != 0.0)[0]
datax = data[nn]

N = datax.shape[0]

N_time = 200

inRes = []
outRes = []

T1 = time.time()

for i in range(0, N, N_time):

  for j in range(i, i+N_time-1):
  
    inRes.append(datax[j, 1:]) # will be the input to the ML model. The order is --> nH, T, nXi
    outRes.append(datax[j+1, 1:]) # will be output of the ML model. The order is --> nH, T, nXi

inRes = np.array(inRes)
outRes = np.array(outRes)

print('Elapsed time = ', time.time() - T1)

print(f'inRes = {inRes}\n\n')
print()
print(f'outRes = {outRes}\n\n')
print(f'inRes.shape = {inRes.shape}')
print()
print(f'outRes.shape = {outRes.shape}')


#---------------------

model = tf.keras.models.load_model('best_model.h5')


new_X = inRes


y_pred = model.predict(new_X)
print("y_pred:", y_pred)

dictx = {'inRes': inRes, 'outRes': outRes, 'y_pred': y_pred}

with open('Test.pkl', 'wb') as f:
  pickle.dump(dictx, f)





