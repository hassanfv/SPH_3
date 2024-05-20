
import pandas as pd
import matplotlib.pyplot as plt

# Load data
file_path = 'CarbonCoolingRates.txt'  # Update with the actual file path
data = pd.read_csv(file_path, delim_whitespace=True, skiprows=19) # Adjust skiprows as necessary

print(data)


# Assign column names based on the information provided
data.columns = ['Temperature', 'C0', 'C+', 'C2+', 'C3+', 'C4+', 'C5+', 'C6+', 'CCIE-only-C']

# Convert data types to float for proper plotting
data = data.astype(float)

# Plotting
plt.figure(figsize=(10, 6))
for col in data.columns[1:]:
    plt.plot(data['Temperature'], data[col], label=col)

# Set logarithmic scale for both axes
plt.xscale('log')
plt.yscale('log')

plt.ylim(-24, -17.4)

plt.xlim(1e4, 1e8)

# Adding labels and title
plt.xlabel('Temperature (K)')
plt.ylabel('Ionic Cooling Efficiency')
plt.title('Ionic Cooling Efficiencies vs Temperature')
plt.legend()

# Show the plot
plt.show()

