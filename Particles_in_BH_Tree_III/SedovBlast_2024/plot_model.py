
import numpy as np
import matplotlib.pyplot as plt

# Read the data from the text files
rho_data = np.loadtxt("rho.dat")
P_data = np.loadtxt("P.dat")
v_data = np.loadtxt("v.dat")
u_data = np.loadtxt("u.dat")
s_data = np.loadtxt("s.dat")

# Separate r_s and the variable values for each dataset
r_s_rho, rho_s = rho_data[:, 0], rho_data[:, 1]
r_s_P, P_s = P_data[:, 0], P_data[:, 1]
r_s_v, v_s = v_data[:, 0], v_data[:, 1]
r_s_u, u_s = u_data[:, 0], u_data[:, 1]
r_s_s, s_s = s_data[:, 0], s_data[:, 1]

# Set up the subplot grid
fig, axs = plt.subplots(2, 3, figsize=(15, 10))

# Plot each variable vs r_s in the respective subplot
axs[0, 0].plot(r_s_rho, rho_s, label='rho_s')
axs[0, 0].set_title('rho_s vs r_s')
axs[0, 0].set_xlabel('r_s')
axs[0, 0].set_ylabel('rho_s')

axs[0, 1].plot(r_s_P, P_s, label='P_s')
axs[0, 1].set_title('P_s vs r_s')
axs[0, 1].set_xlabel('r_s')
axs[0, 1].set_ylabel('P_s')

axs[0, 2].plot(r_s_v, v_s, label='v_s')
axs[0, 2].set_title('v_s vs r_s')
axs[0, 2].set_xlabel('r_s')
axs[0, 2].set_ylabel('v_s')

axs[1, 0].plot(r_s_u, u_s, label='u_s')
axs[1, 0].set_title('u_s vs r_s')
axs[1, 0].set_xlabel('r_s')
axs[1, 0].set_ylabel('u_s')

axs[1, 1].plot(r_s_s, s_s, label='s_s')
axs[1, 1].set_title('s_s vs r_s')
axs[1, 1].set_xlabel('r_s')
axs[1, 1].set_ylabel('s_s')

# Remove the unused subplot (bottom right)
axs[1, 2].axis('off')

# Adjust layout to prevent overlap
plt.tight_layout()

# Display the plot
plt.show()

