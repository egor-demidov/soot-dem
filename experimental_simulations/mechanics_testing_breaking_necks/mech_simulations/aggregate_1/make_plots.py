import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import vtk
from vtk.util.numpy_support import vtk_to_numpy

plt.rcParams.update({'font.size': 14})


strong_path = Path('config_3_strong_focused/stress_strain.csv')
weak_path = Path('config_3_weak/stress_strain.csv')
neck_radius = 8e-9
n_strength = 16.0e9
s_strength = 12.0e9

i1, i2 = 63, 122
files = sorted(
    Path("config_3_strong_focused/run").glob("particles_*.vtk"),
    key=lambda p: int(p.stem.split("_")[-1])
)
positions1 = []
positions2 = []
for file in files:
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(str(file))
    reader.Update()

    points = vtk_to_numpy(reader.GetOutput().GetPoints().GetData())

    positions1.append(points[i1])
    positions2.append(points[i2])

positions1 = np.array(positions1) * 14 # Shape (n_frames, 3)
positions2 = np.array(positions2) * 14 # Shape (n_frames, 3)

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
    f_n[fracture_idx:, :] = 0.0
    f_shear[fracture_idx:, :] = 0.0

    return displacement_mag * 1e9, force_j_mag * 1e9, f_n[:, :] * 1e9, f_shear[:, :] * 1e9

def plot_vectors(ax, v1, v2):
    # Orthonormal basis of the plane (Gram-Schmidt)
    e1 = v1 / np.linalg.norm(v1)

    e2 = v2 - np.dot(v2, e1) * e1
    e2 /= np.linalg.norm(e2)

    # Project onto the plane
    u1 = np.array([np.dot(v1, e1), np.dot(v1, e2)])
    u2 = np.array([np.dot(v2, e1), np.dot(v2, e2)])

    u1 = u1 / np.linalg.norm(u1) * np.log10(np.linalg.norm(u1))
    u2 = u2 / np.linalg.norm(u2) * np.log10(np.linalg.norm(u2))

    print(u1, u2)

    for u, color, label in [
        (u1, 'tab:blue', r'$\mathbf{v}_1$'),
        (u2, 'tab:orange', r'$\mathbf{v}_2$'),
    ]:
        ax.annotate(
            '', xy=u, xytext=(0, 0),
            arrowprops=dict(arrowstyle='->', lw=2, color=color)
        )
    ax.text(*(u * 1.08), label, color=color)

    # Circle at the common origin
    # ax.scatter(0, 0, s=150, facecolors='white',
    #            edgecolors='black', linewidths=2, zorder=5)

    # ax.set_xlim((-4, 4))
    # ax.set_ylim((-4, 4))

    # ax.set_aspect('equal')
    # ax.autoscale()
    # ax.margins(0.2)

disp_strong, force_j_strong, f_n_strong, f_shear_strong = get_force_curve(strong_path)
disp_weak, force_j_weak, f_n_weak, f_shear_weak = get_force_curve(weak_path)

f_n_strong[0, :] = 0.0

# fig, axs = plt.subplot_mosaic([['A', 'B'], ['C', 'D']], figsize=(13, 9.0))
fig, axs = plt.subplot_mosaic([['A', 'B']], figsize=(10, 4.0))


axs['A'].plot(disp_strong, force_j_strong, '-ok')
axs['A'].set_xlim((-1, 5))
axs['A'].set_xlabel('Pulling displacement, nm')
axs['A'].set_ylabel('Pulling force, nN')

axs['B'].plot(disp_weak, force_j_weak, '-ok')
axs['B'].set_xlim((-0.1, 0.5))
axs['B'].set_xlabel('Pulling displacement, nm')
axs['B'].set_ylabel('Pulling force, nN')

fig2, axs2 = plt.subplot_mosaic([['A'], ['B']], figsize=(6.0, 8.0))

axs2['A'].plot(disp_strong, np.linalg.norm(f_shear_strong, axis=1, keepdims=True), '-o', color='tab:orange')
# axs2['A'].axhline(np.pi * neck_radius**2.0 * s_strength * 1e9)
axs2['A'].set_xlim((-1, 5))
axs2['A'].set_xlabel('Pulling displacement, nm')
axs2['A'].set_ylabel('Shear force in neck, nN')

axs2['B'].plot(disp_strong, np.linalg.norm(f_n_strong, axis=1, keepdims=True), '-o', color='tab:blue')
axs2['B'].set_xlim((-1, 5))
# axs2['B'].set_ylim((-0.03, 0.4))
axs2['B'].set_xlabel('Pulling displacement, nm')
axs2['B'].set_ylabel('Normal force in neck, nN')

# plot_vectors(axs['C'], f_n_strong, f_shear_strong)

fig.tight_layout()
fig2.tight_layout()

fig.savefig('aggregate.pdf')
fig2.savefig('neck.pdf')

plt.show()
