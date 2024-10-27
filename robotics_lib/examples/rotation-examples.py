# %% [markdown]
# # Examples of Different Rotation Representations
# This notebook demonstrates various rotation representations and their conversions.

# %%
from sympy import symbols, Matrix, pi, cos, sin, simplify, eye, sqrt, I
from robotics_lib import (
    elem_rot_mat,
    elem_rot_mat_inverse,
    euler_rotation,
    euler_rotation_inverse,
    rpy_rotation,
    angle_axis_to_rot_mat,
    angle_axis_from_rot_mat
)
import numpy as np

# %% [markdown]
# ## 1. Euler Angles (Not Fixed) (e.g., ZX'Z'')
# Example of ZXZ Euler angles sequence

# %%
# Direct computation
# Convert degrees to radians
phi = np.deg2rad(40)
theta = np.deg2rad(60)
psi = np.deg2rad(80)

# Method 1: Using euler_rotation
R_zxz = euler_rotation('zxz', [phi, theta, psi])

# Method 2: Using elementary rotations
R_check = elem_rot_mat('z', phi) @ elem_rot_mat('x', theta) @ elem_rot_mat('z', psi)

# Display results
{
    'euler_rotation_result': R_zxz,
    'elementary_rotations_result': R_check,
    'difference': simplify(R_zxz - R_check)  # Should be zero matrix
}

# %%
# Inverse computation
phi_inv, theta_inv, psi_inv = euler_rotation_inverse('zxz', R_zxz, True)  # True for positive solution

# Convert back to degrees
angles_deg = {
    'phi_deg': float(phi_inv.evalf()) * 180/pi,
    'theta_deg': float(theta_inv.evalf()) * 180/pi,
    'psi_deg': float(psi_inv.evalf()) * 180/pi
}
angles_deg

# %% [markdown]
# ## 2. RPY Angles (Fixed XYZ)
# Example of RPY rotation with fixed axes

# %%
# Direct computation
psi = np.deg2rad(80)    # ROT X
theta = np.deg2rad(60)  # ROT Y
phi = np.deg2rad(40)    # ROT Z

# Method 1: Using rpy_rotation
R_rpy = rpy_rotation('xyz', [psi, theta, phi])

# Method 2: Using elementary rotations (fixed axes)
R_check = elem_rot_mat('z', phi) @ elem_rot_mat('y', theta) @ elem_rot_mat('x', psi)

# Display results
{
    'rpy_rotation_result': R_rpy,
    'elementary_rotations_result': R_check,
    'difference': simplify(R_rpy - R_check)  # Should be zero matrix
}

# %%
# Inverse computation for RPY angles
phi_inv, theta_inv, psi_inv = euler_rotation_inverse('xyz', R_rpy, True)

angles_deg = {
    'phi_deg': float(phi_inv.evalf()) * 180/pi,
    'theta_deg': float(theta_inv.evalf()) * 180/pi,
    'psi_deg': float(psi_inv.evalf()) * 180/pi
}
angles_deg

# %% [markdown]
# ## 3. Angle-Axis Representation
# Example of angle-axis rotation and conversion

# %%
#%% Define axis and angle
r = Matrix([-1, 1, 1])
#%% Normalize the axis
r = r / sqrt(r.dot(r))
#%%
theta = np.deg2rad(65)
theta

#%% Convert to rotation matrix
R = angle_axis_to_rot_mat(r, theta)
R

#%% Convert back to angle-axis
axis_recovered, angle_recovered = angle_axis_from_rot_mat(R)

{
    'original_axis': r,
    'recovered_axis': simplify(axis_recovered),
    'original_angle_deg': 65,
    'recovered_angle_deg': float(angle_recovered.evalf()) * 180/pi
}

# %%
# Singular case example
R_singular = Matrix([
    [-1, 0, 0],
    [0, -1/sqrt(2), -1/sqrt(2)],
    [0, -1/sqrt(2), 1/sqrt(2)]
])

# Try to recover angle-axis
axis_sing, angle_sing = angle_axis_from_rot_mat(R_singular)

{
    'singular_axis': simplify(axis_sing),
    'singular_angle_deg': float(angle_sing.evalf()) * 180/pi
}

# %% [markdown]
# ## 4. Finding Invariant Rotation Axis
# Example of finding the axis that remains unchanged by rotation

# %%
# Create a rotation matrix using euler angles
alpha_1 = -pi/2  # referred to Z
alpha_2 = -pi/4  # referred to Y
alpha_3 = pi/4   # referred to X

#%% R_X(alpha_3)*R_Y(alpha_2)*R_Z(alpha_1)
R_ZYX = euler_rotation('xyz', [alpha_3, alpha_2, alpha_1])

#%% Find eigenvectors and eigenvalues
M = Matrix(R_ZYX)
eigensystem = M.eigenvects()

# Find the eigenvector corresponding to eigenvalue 1
invariant_axis = None
for eigenval, multiplicity, eigenvects in eigensystem:
    if abs(eigenval - 1) < 1e-10:  # Check if eigenvalue is 1
        invariant_axis = eigenvects[0]
        break

# Normalize the invariant axis
if invariant_axis is not None:
    invariant_axis = invariant_axis / sqrt(invariant_axis.dot(invariant_axis))

{
    'rotation_matrix': R_ZYX,
    'invariant_axis': simplify(invariant_axis),
    'verification': simplify(R_ZYX @ invariant_axis - invariant_axis)  # Should be zero vector
}

# %% [markdown]
# ## 5. Combined Example: Multiple Sequence Conversions

# %%
# Start with Euler ZYX sequence
alpha_1 = -pi/2  # Z rotation
alpha_2 = -pi/4  # Y rotation
alpha_3 = pi/4   # X rotation

# Create rotation matrix
R_ZYX = euler_rotation('xyz', [alpha_3, alpha_2, alpha_1])

# Convert to various representations

# 1. Back to Euler angles
euler_angles = euler_rotation_inverse('xyz', R_ZYX, True)

# 2. To angle-axis representation
axis, angle = angle_axis_from_rot_mat(R_ZYX)

# 3. Extract elementary rotations
Rz = elem_rot_mat_inverse('z', R_ZYX)
Ry = elem_rot_mat_inverse('y', elem_rot_mat('z', -Rz) @ R_ZYX)
Rx = elem_rot_mat_inverse('x', elem_rot_mat('y', -Ry) @ elem_rot_mat('z', -Rz) @ R_ZYX)

results = {
    'original_matrix': simplify(R_ZYX),
    'euler_angles': [float(angle.evalf()) for angle in euler_angles],
    'angle_axis': {
        'axis': simplify(axis),
        'angle': float(angle.evalf())
    },
    'elementary_angles': {
        'z_angle': float(Rz.evalf()),
        'y_angle': float(Ry.evalf()),
        'x_angle': float(Rx.evalf())
    }
}
results