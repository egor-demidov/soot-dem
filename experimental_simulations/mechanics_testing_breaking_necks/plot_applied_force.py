from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


force_file = Path('mech_simulations/aggregate_1/config_3/applied_force_log.csv')

data = np.loadtxt(force_file, skiprows=1, delimiter=',')

applied_f_i = data[:, 1:4]
applied_f_j = data[:, 4:7]

net_f_i = data[:, 7:10]
net_f_j = data[:, 10:13]

diff_f_ij = net_f_i - net_f_j

sum_applied_f = applied_f_i + applied_f_j

applied_f_magnitude = np.sqrt(np.sum(applied_f_i * applied_f_i, axis=1))
net_f_i_magnitude = np.sqrt(np.sum(net_f_i * net_f_i, axis=1))
net_f_j_magnitude = np.sqrt(np.sum(net_f_j * net_f_j, axis=1))
diff_f_ij_magnitude = np.sqrt(np.sum(diff_f_ij * diff_f_ij, axis=1))

plt.plot(data[:, 0] * 1e9, applied_f_magnitude * 1e9)
plt.plot(data[:, 0] * 1e9, net_f_i_magnitude * 1e9)
plt.plot(data[:, 0] * 1e9, net_f_j_magnitude * 1e9)
plt.plot(data[:, 0] * 1e9, diff_f_ij_magnitude * 1e9)
plt.xlabel('Time, ns')
plt.ylabel('Force, nN')
plt.show()

