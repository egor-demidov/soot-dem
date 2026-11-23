"""Independent axisymmetric Young–Laplace reference for equal spheres."""
import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares


def bridge_reference(volume, gap, theta_deg, radius=.6, sigma=.008):
    """Solve a symmetric, simple bridge on the facing hemispheres.

    All dimensional inputs must use consistent units. The angle is measured
    through the liquid. An obtuse angle can give an overhang, or a negative
    attractive_force (repulsion). This solver does not establish stability or
    exhaust the set of possible solution branches.
    """
    if not (all(math.isfinite(x) and x > 0 for x in (volume, gap, radius, sigma))
            and 0 < theta_deg < 180):
        raise ValueError("Require positive finite geometry and tension, and 0 < theta < 180")
    volume /= radius**3
    gap /= radius
    theta = math.radians(theta_deg)
    sphere_center = 1 + gap / 2

    def integrate(params, dense=False):
        mid_radius, pressure, length, beta = params

        def rhs(_, state):
            r, z, phi, disk_volume = state
            sin_phi = math.sin(phi)
            return [length * math.cos(phi), length * sin_phi,
                    length * (pressure - sin_phi / max(r, 1e-5)),
                    length * math.pi * r*r * sin_phi]

        def axis_event(_, state):
            return state[0] - 1e-5

        axis_event.terminal = True
        axis_event.direction = -1
        solution = solve_ivp(
            rhs, (0, 1), (mid_radius, 0, math.pi / 2, 0),
            events=axis_event, dense_output=dense,
            rtol=1e-8, atol=1e-10, max_step=.04,
        )
        r, z, phi, disk_volume = solution.y[:, -1]
        cap_height = 1 - math.cos(beta)
        solid_cap = math.pi * cap_height**2 * (1 - cap_height / 3)
        residual = np.array([
            r - math.sin(beta),
            z - gap / 2 - cap_height,
            phi - beta - theta,
            (2 * (disk_volume - solid_cap) - volume) / max(volume, .1),
            1 - solution.t[-1],
        ])
        return residual, solution

    def is_simple(solution):
        r, z, phi, _ = solution.sol(np.linspace(0., 1., 201))
        radial_delta = np.diff(r)
        # A single radial extremum at the midplane admits an axial overhang
        # while excluding looping profiles and interior sphere intersections.
        return (np.min(r) > 0 and np.min(z) >= -1e-7
                and np.min(r*r + (z - sphere_center)**2) >= 1 - 1e-6
                and (np.all(radial_delta >= -1e-8) or np.all(radial_delta <= 1e-8)))

    lower = np.array([1e-4, -500., 1e-4, 1e-4])
    upper = np.array([max(4., 2 * volume**(1/3)), 500.,
                      max(8., 4 * gap), math.pi / 2 - 1e-5])
    best_residual = math.inf
    for beta in (.55, .3, .85, 1.15, .12, 1.45):
        # Start with a circular meridian meeting the sphere at beta+theta.
        # Its signed axial curvature is a useful seed on either side of 90°.
        turn = theta + beta - math.pi / 2
        z = gap / 2 + 1 - math.cos(beta)
        mid = math.sin(beta) + z * math.tan(turn / 2)
        sin_turn = math.sin(turn)
        length = z * (turn / sin_turn if abs(turn) > 1e-6 else 1.)
        pressure = sin_turn / z + 1 / max(mid, .04)
        guess = np.clip([mid, pressure, length, beta], lower + 1e-7, upper - 1e-7)
        fit = least_squares(
            lambda p: integrate(p)[0], guess, bounds=(lower, upper),
            xtol=1e-9, ftol=1e-9, gtol=1e-9, max_nfev=100, x_scale="jac",
        )
        residual = float(np.linalg.norm(fit.fun))
        best_residual = min(best_residual, residual)
        if residual > 1e-6:
            continue
        _, solution = integrate(fit.x, dense=True)
        if not is_simple(solution):
            continue
        mid, pressure, _, beta = map(float, fit.x)
        return dict(
            neck_radius=mid * radius,
            pressure=pressure * sigma / radius,
            attractive_force=math.pi * sigma * radius * (2 * mid - pressure * mid*mid),
            filling_angle_deg=math.degrees(beta),
            residual=residual,
        )
    raise RuntimeError(
        f"No resolved simple bridge on the facing hemispheres at theta={theta_deg:g} deg "
        f"(best residual {best_residual:.3g})"
    )
