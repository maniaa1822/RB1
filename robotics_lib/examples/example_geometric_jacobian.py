#%% Import required libraries
from sympy import symbols, pi, Matrix, simplify
from robotics_lib.kinematics import dh_matrix, get_f_r, geometric_jacobian
import numpy as np
import matplotlib.pyplot as plt
#%matplotlib inline

#%% Define symbolic variables
q1, q2, q3 = symbols('q1 q2 q3', real=True)
l1 = symbols('l1', real=True)

#%% Create DH parameter table for a RRR robot
# Example: 3-DOF robot with all revolute joints
dh_table = Matrix([
    [0,      l1,     0,  q1],           # alpha, a, d, theta for joint 1
    [pi/2,   0,      0,  q2 + pi/2],    # alpha, a, d, theta for joint 2
    [0,      0,      q3, 0]             # alpha, a, d, theta for joint 3
])

dh_table

#%% Calculate transformation matrices and forward kinematics
T, A = dh_matrix(dh_table)
f_r_3D = get_f_r(T)
simplify(f_r_3D)

#%% Calculate geometric Jacobian
# Define joint sequence ('r' for revolute, 'p' for prismatic)
sequence = 'rrp'  # our robot has two revolute joints and one prismatic

# Define joint variables
q_vars = Matrix([q1, q2, q3])

#%% Calculate geometric Jacobian
Jl, Ja = geometric_jacobian(f_r_3D, sequence, q_vars, dh_table)

#%% Display linear part
print("Linear part of geometric Jacobian:")
simplify(Jl)

#%% Display angular part
print("\nAngular part of geometric Jacobian:")
simplify(Ja)

#%% Verify Jacobian at specific configuration
# Let's substitute some values
config = {q1: pi/4, q2: pi/3, q3: 0.5, l1: 1.0}

# Calculate Jacobian at this configuration
Jl_num = Jl.subs(config)
Ja_num = Ja.subs(config)

# Combine into full geometric Jacobian
J_geometric = Matrix.vstack(Jl_num, Ja_num)

# Display numerical result
J_geometric.evalf()
