#%% Import required libraries
from sympy import symbols, pi, Matrix, simplify
from robotics_lib.kinematics import dh_matrix, get_f_r, get_rotation_mat
import numpy as np
import matplotlib.pyplot as plt
#%matplotlib inline
#%% Define symbolic variables
q1, q2, q3 = symbols('q1 q2 q3', real=True)
l1 = symbols('l1', real=True)

#%% Create DH parameter table
dh_table = Matrix([
    [0,      l1,     0,  q1],           # alpha, a, d, theta for joint 1
    [pi/2,   0,      0,  q2 + pi/2],    # alpha, a, d, theta for joint 2
    [0,      0,      q3, 0]             # alpha, a, d, theta for joint 3
])

# Display DH table
dh_table

#%% Calculate transformation matrices
# Get final transformation matrix T and individual matrices A
T, A = dh_matrix(dh_table)

# Display final transformation matrix
simplify(T)

#%% Display individual transformation matrices
# A0_1
simplify(A[0])

# A1_2
simplify(A[1])

# A2_3
simplify(A[2])

#%% Calculate forward kinematics
# Get the forward kinematics mapping [X Y Z PHI]
f_r_3D = get_f_r(T)
simplify(f_r_3D)

# Select only [X Y PHI] components
f_r = Matrix([f_r_3D[0], f_r_3D[1], f_r_3D[3]])
simplify(f_r)

#%% Calculate Jacobian
vars_list = Matrix([q1, q2, q3])
j = f_r.jacobian(vars_list)
simplify(j)

#%% Transform Jacobian to frame 1
R0_1 = get_rotation_mat(A[0])
j1_3 = R0_1.transpose() * j
simplify(j1_3)

#%% Calculate determinant
det_j1_3 = simplify(j1_3.det())
det_j1_3

#%% Calculate torques for specific configurations
# Define force vector
f = Matrix([0, 1.5, -4.5])

# Configuration 1: q1 = π/2, q2 = 0, q3 = 3, l1 = 0.5
j_0 = j.subs({q1: pi/2, q2: 0, q3: 3, l1: 0.5})
tau_0 = -j_0 * f
tau_0.evalf()

# Configuration 2: q1 = 0, q2 = π/2, q3 = 0, l1 = 0.5
j_s = j.subs({q1: 0, q2: pi/2, q3: 0, l1: 0.5})
tau_s = -j_s * f
tau_s.evalf()

#%% Visualize results
plt.style.use('seaborn')
plt.figure(figsize=(12, 8))

# Calculate torques for different q1 values
q1_vals = np.linspace(-pi/2, pi/2, 100)
torques = []

for q1_val in q1_vals:
    j_val = j.subs({q1: q1_val, q2: 0, q3: 0, l1: 0.5})
    tau = -j_val * f
    torques.append([float(tau[i].evalf()) for i in range(3)])

torques = np.array(torques)

# Create plot
plt.plot(q1_vals, torques[:, 0], label='τ1', linewidth=2)
plt.plot(q1_vals, torques[:, 1], label='τ2', linewidth=2)
plt.plot(q1_vals, torques[:, 2], label='τ3', linewidth=2)
plt.grid(True, alpha=0.3)
plt.xlabel('Joint angle q1 (rad)', fontsize=12)
plt.ylabel('Torque (N⋅m)', fontsize=12)
plt.title('Joint Torques vs q1', fontsize=14, pad=20)
plt.legend(fontsize=10)
plt.tight_layout()