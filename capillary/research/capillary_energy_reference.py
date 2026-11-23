"""Independent geometric oracle for capillary virtual work; NOT an SPH solver.

Area energies on a supplied triangle mesh make force/torque identities testable
without assuming an interface reconstruction or a fluid/solid force partition.
The production adapter does not import this module.
"""
import numpy as np


def area_energy_forces(points, triangles, coefficients):
    """Return E=sum gamma_t A_t and all -dE/dx, with analytic area derivatives.

    A negative coefficient is allowed for the solid-liquid excess energy.
    Every vertex dependence is differentiated; no force balancing is applied.
    """
    points = np.asarray(points, dtype=float)
    triangles = np.asarray(triangles, dtype=int)
    coefficients = np.broadcast_to(np.asarray(coefficients, dtype=float), (len(triangles),))
    if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
        raise ValueError('Expected finite 3D vertices')
    if triangles.ndim != 2 or triangles.shape[1] != 3 or not np.isfinite(coefficients).all():
        raise ValueError('Expected triangles and finite energy coefficients')
    if np.any(triangles < 0) or np.any(triangles >= len(points)):
        raise ValueError('Triangle index out of range')
    a,b,c = (points[triangles[:,i]] for i in range(3))
    cross = np.cross(b-a,c-a)
    norm = np.linalg.norm(cross,axis=1)
    if np.any(norm == 0): raise ValueError('Degenerate interface triangle')
    normal = cross/norm[:,None]
    forces = np.zeros_like(points)
    for index,edge in enumerate((b-c,c-a,a-b)):
        gradient = .5*np.cross(edge,normal)
        np.add.at(forces,triangles[:,index],-coefficients[:,None]*gradient)
    return float(np.sum(.5*coefficients*norm)),forces


def rigid_load(points, forces, attached, center):
    """Force and torque on one rigid body owning the selected interface vertices."""
    relative = np.asarray(points)[attached]-np.asarray(center)
    selected = np.asarray(forces)[attached]
    return selected.sum(axis=0),np.cross(relative,selected).sum(axis=0)


def spherical_cap(radius, theta_deg, latitude_rings=24, azimuthal_points=96):
    """Cap plus wetted disk. Shared contact vertices have z=0.

    Contact vertices are marked as solid-attached only to evaluate a frozen
    contact geometry's rigid variation. Real wetting also needs contact-line
    slip degrees of freedom; permanent attachment would pin it artificially.
    """
    if not (radius>0 and 0<theta_deg<180 and latitude_rings>=2 and azimuthal_points>=3):
        raise ValueError('Invalid cap geometry')
    theta=np.deg2rad(theta_deg)
    zc=-radius*np.cos(theta)
    points=[[0.,0.,zc+radius]]
    for phi in np.linspace(theta/latitude_rings,theta,latitude_rings):
        for angle in np.arange(azimuthal_points)*2*np.pi/azimuthal_points:
            points.append([radius*np.sin(phi)*np.cos(angle),radius*np.sin(phi)*np.sin(angle),zc+radius*np.cos(phi)])
    top=[]
    for k in range(azimuthal_points): top.append([0,1+k,1+(k+1)%azimuthal_points])
    for ring in range(latitude_rings-1):
        first=1+ring*azimuthal_points;second=first+azimuthal_points
        for k in range(azimuthal_points):
            nxt=(k+1)%azimuthal_points
            top.extend([[first+k,second+k,second+nxt],[first+k,second+nxt,first+nxt]])
    boundary=np.arange(1+(latitude_rings-1)*azimuthal_points,1+latitude_rings*azimuthal_points)
    disk_center=len(points);points.append([0.,0.,0.])
    disk=[[disk_center,int(boundary[(k+1)%azimuthal_points]),int(boundary[k])] for k in range(azimuthal_points)]
    attached=np.r_[boundary,disk_center]
    return np.asarray(points),np.asarray(top),np.asarray(disk),attached
