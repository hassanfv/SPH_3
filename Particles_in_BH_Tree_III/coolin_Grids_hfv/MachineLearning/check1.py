
import numpy as np
import pickle
import pandas as pd
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt


with open('chimesResLOG.pkl', 'rb') as f:
  data = pickle.load(f)

#print(data.keys())
#['t_Arr_in_yrs', 'TEvol', 'nHe0', 'nHep', 'nHepp', 'nH0', 'nHp', 'nC0', 'nC1', 'nC2', 'nC3', 'nC4', 'nC5', 'nC6']

t_yrs = data['t_Arr_in_yrs']
T = data['TEvol']

nH0 = data['nH0']
nH0 = data['nHp']

nH0 = data['nHe0']
nH0 = data['nHep']
nH0 = data['nHepp']

nC0 = data['nC0']
nC1 = data['nC1']
nC2 = data['nC2']
nC3 = data['nC3']
nC4 = data['nC4']
nC5 = data['nC5']
nC6 = data['nC6']

print(t_yrs.shape, T.shape, nH0.shape, nC2.shape)

df = pd.DataFrame(data)

df = df.drop('t_Arr_in_yrs', axis = 1)
print(df.head())

X = df[:-1]  # all rows except the last one
y = df[1:]   # all rows except the first one

X_train = X_test = X #, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
y_train = y_test = y



# Scale the data
#scaler = StandardScaler()
#X_train = scaler.fit_transform(X_train)
#X_test = scaler.transform(X_test)


# Build the model
model = Sequential([
    Dense(500, activation='relu', input_dim=X_train.shape[1]),
    Dense(200, activation='relu'),
    Dense(100, activation='relu'),
    Dense(100, activation='relu'),
    Dense(y_train.shape[1], activation='linear')  # or 'sigmoid', depending on your specific case
])

# Compile the model
model.compile(optimizer='adam', loss='mean_squared_error')

# Train the model
model.fit(X_train, y_train, epochs=10, batch_size=32)

# Evaluate the model
loss = model.evaluate(X_test, y_test)
print('Test loss:', loss)

y_pred = model.predict(X_test)

# Reverse the scaling of predictions
#y_pred_orig = scaler.inverse_transform(y_pred)
#y_test_orig = scaler.inverse_transform(y_test)

print(y_pred.shape, y_test.shape)

plt.scatter(y_test.iloc[:, 0], y_pred[:, 0], s = 1, color = 'k')
plt.plot(y_test.iloc[:, 0], y_test.iloc[:, 0], color = 'r')

plt.show()



