
import pandas as pd
import matplotlib.pyplot as plt

file_path = 'CarbonCoolingRates.txt'
df = pd.read_csv(file_path, delim_whitespace=True, skiprows=19) # Adjust skiprows as necessary
df.columns = ['Temp', 'C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'CCIE-only-C']

df.to_csv('CData.csv', index = False)

Temp = df['Temp'].values
C0 = df['C0'].values
C1 = df['C1'].values
C2 = df['C2'].values
C3 = df['C3'].values
C4 = df['C4'].values
C5 = df['C5'].values
C6 = df['C6'].values

# Plotting
plt.figure(figsize=(10, 6))

plt.plot(Temp, C0, label = 'C0')
plt.plot(Temp, C1, label = 'C1')
plt.plot(Temp, C2, label = 'C2')
plt.plot(Temp, C3, label = 'C3')
plt.plot(Temp, C4, label = 'C4')
plt.plot(Temp, C5, label = 'C5')
plt.plot(Temp, C6, label = 'C6')

# Set logarithmic scale for both axes
plt.xscale('log')
plt.yscale('log')

plt.ylim(-24, -17.4)
plt.xlim(9e3, 2e8)

# Adding labels and title
plt.xlabel('Temperature (K)')
plt.ylabel('Ionic Cooling Efficiency')
plt.title('Ionic Cooling Efficiencies vs Temperature')
plt.legend()

# Show the plot
plt.show()

