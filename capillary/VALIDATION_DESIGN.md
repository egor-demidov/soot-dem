# Reference equations and diagnostic interpretation

## Droplet
A spherical three-dimensional liquid interface has delta-p = 2 sigma/R.
For a spherical cap on z=0 with sphere center z_c, measured through liquid,
cos(theta) = -z_c/R. The measurements fit positions, never the interface normals
modified by the wetting routine. A finite particle spacing displaces the sampled
interface; raw and half-spacing-offset estimates are both reported.

## Equal-sphere Young–Laplace reference
Two sphere centers are at z=+/-(R+gap/2). The equilibrium bridge is symmetric
about z=0. Follow its meridian from the middle to the upper contact line by arc
length s. Let phi be the tangent angle to the radial axis:

    dr/ds = cos(phi)
    dz/ds = sin(phi)
    dphi/ds = delta-p/sigma - sin(phi)/r
    dV_half_gross/ds = pi r^2 sin(phi)

At the middle: r=r_neck, z=0, phi=pi/2, V_half_gross=0.
At the upper contact line of filling angle beta:

    r = R sin(beta)
    z = gap/2 + R (1-cos(beta))
    phi = beta + theta
    V = 2 [V_half_gross - pi h^2 (R-h/3)], h=R(1-cos(beta)).

The code solves for r_neck, delta-p, half arc length and beta simultaneously.
The tensile axial force (positive attraction magnitude) is

    F = 2 pi sigma r_neck - pi delta-p r_neck^2.

This is a separate boundary-value calculation, not the force used to advance
solids in the SPH case. An exact cylinder, beta+theta=pi/2, verifies pressure
sigma/r and force pi sigma r. Similarity under length scaling is also tested.
Multiple branches can exist; root convergence alone does not establish branch
stability. The reference solver is restricted to 0<theta<90 degrees and a
resolved symmetric bridge. Final SPH shape and quasi-static state must be checked
before making a force comparison. An offset initial cylinder tests relaxation;
its transient asymmetric force is not compared to this equilibrium law.

## Supported mechanical balance
For a fixed outer container and any fixed spheres, the fluid plus free spheres
is an open mechanical subsystem. Let I_support be the integrated force received
by all fixed solids and J_support its moment about the fixed origin. Test

    P_fluid + P_free_spheres + I_support = initial value
    L_fluid + L_free_spheres + J_support = initial value.

For dry spheres: P=mV and L=X cross mV + I omega, I=2mR^2/5. A spherical inertia
is independent of orientation, so no omitted body-frame gyroscopic term occurs.
Quaternion integration rotates surface samples consistently. Forward impulse and
position splitting are first order; the half-step test separates this error from
spatial and force-transfer errors. No correction is applied to make these
balances pass. Pressure extrapolation includes fluid acceleration, so an audit
must not assume all surface-tension force is absent from the wall reaction.

The moving case is intentionally a coupling diagnostic and stops at contact. It
is not the final DEM integrator or an aggregate simulation. Its equilibrium bridge
comparison, dynamic force validation and torque accuracy are separate questions.
