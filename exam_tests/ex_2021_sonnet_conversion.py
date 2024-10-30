#%% Imports and Setup 
import numpy as np
import sympy as sp
from sympy import Matrix, cos, sin, simplify, atan2, pi, symbols
import matplotlib.pyplot as plt
from math import radians, degrees

# Import our robotics library functions
from robotics_lib import *
# Define symbolic variables
q1, q2, q3, q4 = symbols('q1 q2 q3 q4', real=True)
l1, l2, l3 = symbols('l1 l2 l3', real=True)
alpha, beta, gamma = symbols('alpha beta gamma', real=True)
t, T = symbols('t T', real=True)

def S(v):
    """Create skew-symmetric matrix"""
    return Matrix([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])

#%% YXZ Euler Rotation
R_yxz = euler_rotation('yxz', [alpha, beta, gamma])
R_yxz

#%% Calculate angles
b = float(atan2(-np.sqrt(3)/2, 0.5))
a = float(atan2(0, -0.5/cos(b)))
c = float(atan2(0.5/cos(b), 0))

(f"a = {degrees(a):.4f}°", 
 f"b = {degrees(b):.4f}°", 
 f"c = {degrees(c):.4f}°")

#%% Substitute values and calculate WR_B
R = R_yxz.subs({
    alpha: radians(180),
    beta: radians(-60),
    gamma: radians(90)
})
WR_B = R

# Display both matrices
('R =', R, 'WR_B =', WR_B)

#%% Angular velocity
BV_B = Matrix([1, -1, 0])
WR_B_dot = S(WR_B * BV_B) * WR_B

('BV_B =', BV_B, 'WR_B_dot =', WR_B_dot)

#%% First Robot Configuration - DH Parameters
alpha = [-pi/2, -pi/2, 0]
a = [-l1, -l2, -l3]
d = [symbols('l0'), 0, 0]
theta = [q1, q2, q3]

table = Matrix([
    [alpha[0], a[0], d[0], theta[0]],
    [alpha[1], a[1], d[1], theta[1]],
    [alpha[2], a[2], d[2], theta[2]]
])

('DH Parameters:', table)

#%% Calculate transformations
T, A = dh_matrix(table)

('A0_1 =', A[0], 
 'A1_2 =', A[1], 
 'A2_3 =', A[2], 
 'T =', T)

#%% Forward kinematics
f_r_3D = get_f_r(T)
f_r = Matrix([f_r_3D[0], f_r_3D[1], f_r_3D[2]])

('f_r_3D =', f_r_3D, 'f_r =', f_r)

#%% Second Robot Configuration
alpha2 = [0, -pi/2, pi/2, 0]
a2 = [l1, 0, 0, l2]
d2 = [0, 0, q3, 0]
theta2 = [q1, q2-pi/2, 0, q4+pi/2]

table2 = Matrix([
    [alpha2[0], a2[0], d2[0], theta2[0]],
    [alpha2[1], a2[1], d2[1], theta2[1]],
    [alpha2[2], a2[2], d2[2], theta2[2]],
    [alpha2[3], a2[3], d2[3], theta2[3]]
])

('Second Robot DH Parameters:', table2)

#%% Second Robot Forward Kinematics and Jacobian
T2, A2 = dh_matrix(table2)
f_r_3D2 = get_f_r(T2)
f_r2 = Matrix([f_r_3D2[0], f_r_3D2[1], f_r_3D2[3]])  # Using alpha_z instead of z
j = simplify(f_r2.jacobian([q1, q2, q3, q4]))

('T2 =', T2,
 'f_r2 =', f_r2,
 'Jacobian =', j)

#%% Planar 2R Robot
l1_val = 1.0
l2_val = 1.0

def planar_2R_kinematics(q1, q2, l1=l1_val, l2=l2_val):
    r = Matrix([
        l1*cos(q1) + l2*cos(q1+q2),
        l1*sin(q1) + l2*sin(q1+q2)
    ])
    J = r.jacobian([q1, q2])
    return r, J

r, J = planar_2R_kinematics(q1, q2)
('r =', r, 'J =', J)

#%% Test configurations
qa = Matrix([3*pi/4, -pi/2])
qb = Matrix([pi/2, -pi/2])
FA = (10/np.sqrt(2))*Matrix([1, 1])

J_qa = J.subs({q1: qa[0], q2: qa[1]})
tau_a = J_qa.transpose() * FA

FB = (10/np.sqrt(2))*Matrix([-1, -1])
J_qb = J.subs({q1: qb[0], q2: qb[1]})
tau_b = -J_qb.transpose() * FB

('tau_a =', tau_a, 'tau_b =', tau_b)

#%% Trajectory Generation
def cubic_polynomial(t, T, qi, qf, vi, vf):
    """Generate cubic polynomial trajectory"""
    Dq = qf - qi
    
    # Calculate coefficients
    a = (T*vf + T*vi - 2*Dq)/Dq
    b = (3*Dq - 2*T*vi - T*vf)/Dq
    c = T*vi/Dq
    d = 0
    
    # Normalized time
    tau = t/T
    
    # Generate trajectory
    q = qi + Dq*(a*tau**3 + b*tau**2 + c*tau + d)
    return q

#%% Joint 1 trajectory coefficients
T_val = 2
qi = pi
qf = 1.4805
vi = 0
vf = -0.7185
Dq = qf - qi

a1 = (T_val*vf + T_val*vi - 2*Dq)/Dq
b1 = (3*Dq - 2*T_val*vi - T_val*vf)/Dq
c1 = T_val*vi/Dq
d1 = 0

('Joint 1 coefficients:',
 f'a = {a1:.4f}',
 f'b = {b1:.4f}',
 f'c = {c1:.4f}',
 f'd = {d1:.4f}')

#%% Joint 2 coefficients
qi = 0
qf = -2.3016
vi = 0
vf = 0.4601
Dq = qf - qi

a2 = (T_val*vf + T_val*vi - 2*Dq)/Dq
b2 = (3*Dq - 2*T_val*vi - T_val*vf)/Dq
c2 = T_val*vi/Dq
d2 = 0

('Joint 2 coefficients:',
 f'a = {a2:.4f}',
 f'b = {b2:.4f}',
 f'c = {c2:.4f}',
 f'd = {d2:.4f}')

#%% Generate and plot trajectories
t_vals = np.linspace(0, 2, 100)
q_traj1 = [cubic_polynomial(t, 2, pi, 1.4805, 0, -0.7185) for t in t_vals]
q_traj2 = [cubic_polynomial(t, 2, 0, -2.3016, 0, 0.4601) for t in t_vals]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

ax1.plot(t_vals, q_traj1)
ax1.set_title('Joint 1 Trajectory')
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Position (rad)')
ax1.grid(True)

ax2.plot(t_vals, q_traj2)
ax2.set_title('Joint 2 Trajectory')
ax2.set_xlabel('Time (s)')
ax2.set_ylabel('Position (rad)')
ax2.grid(True)

plt.tight_layout()
plt.show()
# %%
