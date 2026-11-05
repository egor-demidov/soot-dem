import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path


strong_path = Path('config_3_strong_focused/stress_strain.csv')
weak_path = Path('config_3_weak/stress_strain.csv')

def get_force_curve(data_path):
    data = np.loadtxt(data_path, delimiter=',', skiprows=1)

    displacement = data[:, 1:4] - data[0, 1:4]
    force_j = data[:, 7:10]
    f_n = data[:, 10:13]
    f_t = data[:, 13:16]
    f_r = data[:, 16:19]
    f_o = data[:, 19:]

    f_shear = f_t + f_r + f_o

    displacement_mag = np.linalg.norm(displacement, axis=1, keepdims=True)
    force_j_mag = np.linalg.norm(force_j, axis=1, keepdims=True)
    fracture_idx = np.argmax(force_j_mag)

    print(f_n[fracture_idx, :])
    print(f_shear[fracture_idx, :])

    return displacement_mag, force_j_mag


disp_strong, force_j_strong = get_force_curve(strong_path)
disp_weak, force_j_weak = get_force_curve(weak_path)

fig, axs = plt.subplot_mosaic([['A', 'B']], figsize=(13, 4.0))


axs['A'].plot(disp_strong, force_j_strong, '-o')
axs['A'].set_xlim((-1e-9, 1e-8))
axs['A'].set_xlabel('Displacement')
axs['A'].set_ylabel('Force')

axs['B'].plot(disp_weak, force_j_weak, '-o')
axs['B'].set_xlim((-1e-10, 1e-9))
axs['B'].set_xlabel('Displacement')
axs['B'].set_ylabel('Force')

plt.show()
