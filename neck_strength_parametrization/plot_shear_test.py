import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


INP_PATH = Path('../cmake-build-debug/shear_forces.csv')

COL_LABELS = ['n', 'dx', 'dy', 'dz', 'xi_n', 'xi_tx', 'xi_ty', 'xi_tz', 'xi_rx', 'xi_ry', 'xi_rz', 'xi_ox', 'xi_oy', 'xi_oz', 'fx', 'fy', 'fz', 'taux', 'tauy', 'tauz']
COL_MAP = {k: i for i, k in enumerate(COL_LABELS)}

data = np.loadtxt(INP_PATH, delimiter=',', skiprows=1)

fig, ax = plt.subplots()
ax2 = ax.twinx()

# ax.plot(data[:, COL_MAP['dx']], data[:, COL_MAP['fx']])

ax2.plot(data[:, COL_MAP['dx']], data[:, COL_MAP['xi_n']], label='Xn')
ax2.plot(data[:, COL_MAP['dx']], data[:, COL_MAP['xi_tx']], label='Xtx')
ax2.plot(data[:, COL_MAP['dx']], data[:, COL_MAP['xi_tz']], label='Xtz')

ax2.legend()

plt.show()
