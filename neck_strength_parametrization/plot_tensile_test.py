import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


INP_PATH = Path('../cmake-build-debug/tensile_forces.csv')

COL_LABELS = ['n', 'dx', 'dy', 'dz', 'xi_n', 'xi_tx', 'xi_ty', 'xi_tz', 'xi_rx', 'xi_ry', 'xi_rz', 'xi_ox', 'xi_oy', 'xi_oz', 'fx', 'fy', 'fz', 'taux', 'tauy', 'tauz']
COL_MAP = {k: i for i, k in enumerate(COL_LABELS)}

data = np.loadtxt(INP_PATH, delimiter=',', skiprows=1)

plt.plot(data[:, COL_MAP['dz']] - data[0, COL_MAP['dz']], data[:, COL_MAP['fz']])
plt.show()
