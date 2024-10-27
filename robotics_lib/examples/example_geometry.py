# %% [markdown]
# # Complete Examples for Geometric Functions Library
# This notebook provides examples for all functions in the geometric library.

# %%
from sympy import symbols, Matrix, pi, cos, sin, simplify, eye, sqrt
from robotics_lib.geometric import *

# %% [markdown]
# ## 1. Elementary Rotation Matrices (elem_rot_mat)
# Examples of rotation matrices around each axis, both symbolic and numerical.

# %%
# Create symbolic angles
theta, alpha, beta = symbols('theta alpha beta')

# Rotations around each axis
Rx = elem_rot_mat('x', alpha)
Ry = elem_rot_mat('y', beta)
Rz = elem_rot_mat('z', theta)

# Display symbolic results
{
    'Rx': Rx,
    'Ry': Ry,
    'Rz': Rz
}

# %%
# Numerical example: evaluate for specific angles
angles_rad = {alpha: pi/4, beta: pi/3, theta: pi/6}
{
    'Rx_num': Rx.subs(angles_rad).evalf(),
    'Ry_num': Ry.subs(angles_rad).evalf(),
    'Rz_num': Rz.subs(angles_rad).evalf()
}

# %% [markdown]
# ## 2. Euler Rotation Matrices (euler_rotation)
# Examples of different Euler angle sequences.

# %%
# Define symbolic Euler angles
phi, theta, psi = symbols('phi theta psi')

# Different Euler sequences
R_xyz = euler_rotation('xyz', [phi, theta, psi])
R_zyz = euler_rotation('zyz', [phi, theta, psi])
R_zyx = euler_rotation('zyx', [phi, theta, psi])

# Display symbolic results
{
    'xyz_sequence': simplify(R_xyz),
    'zyz_sequence': simplify(R_zyz),
    'zyx_sequence': simplify(R_zyx)
}

# %%
# Numerical example
euler_angles = {phi: pi/6, theta: pi/4, psi: pi/3}
R_xyz.subs(euler_angles).evalf()

# %% [markdown]
# ## 3. Roll-Pitch-Yaw Rotations (rpy_rotation)
# Examples of RPY rotations with different sequences.

# %%
# Define symbolic RPY angles
roll, pitch, yaw = symbols('roll pitch yaw')

# Different RPY sequences
R_rpy_xyz = rpy_rotation('xyz', [roll, pitch, yaw])
R_rpy_zyx = rpy_rotation('zyx', [roll, pitch, yaw])

# Display symbolic results
{
    'rpy_xyz': simplify(R_rpy_xyz),
    'rpy_zyx': simplify(R_rpy_zyx)
}

# %%
# Compare RPY with equivalent Euler rotation
# RPY xyz is equivalent to Euler zyx with reversed angles
diff = simplify(R_rpy_xyz - euler_rotation('zyx', [yaw, pitch, roll]))
diff  # Should be zero matrix

# %% [markdown]
# ## 4. Elementary Rotation Matrix Inverse
# Example of extracting rotation angle from rotation matrix.

# %%
# Create a symbolic rotation matrix
theta = symbols('theta')
R = elem_rot_mat('z', theta)
#%%
# Extract the angle back
theta_recovered = elem_rot_mat_inverse('z', R)
# Should return theta
simplify(theta_recovered)

# %% [markdown]
# ## 5. Euler Rotation Matrix Inverse
# Example of extracting Euler angles from rotation matrix.

# %%
# Create a symbolic rotation matrix using Euler angles
phi, theta, psi = symbols('phi theta psi')
R_euler = euler_rotation('zyx', [phi, theta, psi])

# Extract angles back
phi_rec, theta_rec, psi_rec = euler_rotation_inverse('zyx', R_euler)
{
    'phi': simplify(phi_rec),
    'theta': simplify(theta_rec),
    'psi': simplify(psi_rec)
}

# %% [markdown]
# ## 6. Angle-Axis Representations
# Example of converting between angle-axis and rotation matrix representations.

# %%
# Define symbolic axis and angle
rx, ry, rz, alpha = symbols('rx ry rz alpha')
r = Matrix([rx, ry, rz])

# Convert from angle-axis to rotation matrix
R_aa = angle_axis_to_rot_mat(r, alpha)
R_aa

# %%
# Convert back to angle-axis representation
axis_recovered, angle_recovered = angle_axis_from_rot_mat(R_aa)
{
    'recovered_axis': simplify(axis_recovered),
    'recovered_angle': simplify(angle_recovered)
}

# %%
# Numerical example: rotation of pi/3 around [1,1,1]
subs_dict = {
    rx: 1/sqrt(3),  # Normalized vector
    ry: 1/sqrt(3),
    rz: 1/sqrt(3),
    alpha: pi/3
}
R_aa.subs(subs_dict).evalf()

# %% [markdown]
# ## 5. Homogeneous Transformations (affine_T, affine_get_translation, affine_get_R)
# Examples of creating and manipulating homogeneous transformation matrices.

# %%
# Create symbolic variables
theta = symbols('theta')
tx, ty, tz = symbols('tx ty tz')

# Create rotation and translation
R = elem_rot_mat('z', theta)
t = Matrix([tx, ty, tz])

# Create homogeneous transformation matrix
T = affine_T(R, t)
T

# %%
# Extract components
{
    'rotation_matrix': affine_get_R(T),
    'translation_vector': affine_get_translation(T)
}

# %%
# Numerical example
subs_dict = {theta: pi/4, tx: 1, ty: 2, tz: 3}
T_num = T.subs(subs_dict).evalf()
{
    'full_transform': T_num,
    'extracted_rotation': affine_get_R(T_num).evalf(),
    'extracted_translation': affine_get_translation(T_num).evalf()
}
