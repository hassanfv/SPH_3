
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from tensorflow.keras.callbacks import ModelCheckpoint
from sklearn.model_selection import train_test_split

with open ("initialData.pkl", "rb") as f:
  data = pickle.load(f)

#print(data.shape)

nx = np.where(data[:, 0] == 0.0)
#print(nx)

nn = np.where(data[:, 0] != 0.0)[0]
datax = data[nn]


#print(datax.shape)
#print()
#print(datax)

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

print()
print(np.sort(inRes[:, 1]))

print(f'inRes = {inRes}\n\n')
print()
print(f'outRes = {outRes}\n\n')
print(f'inRes.shape = {inRes.shape}')
print()
print(f'outRes.shape = {outRes.shape}')








X = inRes
Y = outRes


print(X[0])


M = Y.shape[1]

X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2, random_state=42)

model = Sequential([
    Dense(256, activation='relu', input_shape=(X_train.shape[1],)),  # First hidden layer
    Dense(128, activation='relu'),  # Second hidden layer
    Dense(64, activation='relu'),  # Additional third hidden layer
    Dense(M)  # Output layer with M neurons for regression
])


# Compile the model
model.compile(optimizer='adam', loss='mean_squared_error')

# Setup ModelCheckpoint
checkpoint_filepath = 'best_model.h5'
model_checkpoint_callback = ModelCheckpoint(
    filepath=checkpoint_filepath,
    save_best_only=True,
    monitor='val_loss',
    mode='min',
    verbose=1
)

# Train the model
history = model.fit(
    X_train, Y_train,
    epochs=5,
    validation_split=0.2,
    callbacks=[model_checkpoint_callback],
    verbose=1
)

# You can also evaluate to print the test loss after training if needed
test_loss = model.evaluate(X_test, Y_test)
print(f'Test Loss: {test_loss}')






