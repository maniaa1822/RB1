# rotation with euler angles ex 2021

#%%import 
from robotics_lib.geometric import *
from robotics_lib.kinematics import *
from sympy import symbols, Matrix, pi, cos, sin, simplify, eye, sqrt,deg,solve
from robotics_lib.examples.rotation-examples import R_singular


#%%
alpha, beta, theta = symbols('alpha beta theta')

#%% euler_rotation
R_yxz = euler_rotation('yxz', [alpha, beta, theta])
R_yxz
# %% compute b, a, c with atan2
b = atan2(-sqrt(3)/2, 1/2)
a = atan2(0, (1/2)/cos(b))
c = atan2((1/2)/cos(b),0)


# %% convert to degrees the output should be a = 180, b = -60 , c = 90
a = deg(a)
b = deg(b)
c = deg(c)

#%% get numerical values
a_num = a.evalf()
b_num = b.evalf()
c_num = c.evalf()

# %% display results
{
    'a': a,
    'b': b,
    'c': c
}


# %% sub in the values of a, b, c to get the rotation matrix
R_yxz_num = R_yxz.subs({alpha: a, beta: b, theta: c})

#%% create a new vector BV_B = [1, -1, 0]
BV_B = Matrix([1, -1, 0])
BV_B

WR_B = R_yxz_num
WR_B

#%%WR_B * BV_B 
WR_B_BV_B = WR_B @ BV_B
WR_B_BV_B

#%% calculate the skew symmetric matrix of WR_B_BV_B
K = Matrix([
    [0, -WR_B_BV_B[2], WR_B_BV_B[1]],
    [WR_B_BV_B[2], 0, -WR_B_BV_B[0]],
    [-WR_B_BV_B[1], WR_B_BV_B[0], 0]
])

#%% multiply WR_B by K
WR_B_dot = K @ WR_B
WR_B_dot

# second part of the question

#create q1, q2, q3 symbols and l0, l1, l2, l3 symbols
q1, q2, q3 = symbols('q1 q2 q3')
l0, l1, l2, l3 = symbols('l0 l1 l2 l3')

# %% create the DH table
dh_table = Matrix([
    [-pi/2, -l1, l0, q1],           # alpha, a, d, theta for joint 1
    [-pi/2, -l1, 0, q2],    # alpha, a, d, theta for joint 2
    [0, -l3, 0, q3]             # alpha, a, d, theta for joint 3
])

dh_table
# %%
T, A = dh_matrix(dh_table)
# %% display final transformation matrix and individual transformation matrices assinging
#individual transformation matrices to A0_1, A1_2, A2_3
T = simplify(T)
A0_1 = simplify(A[0])
A1_2 = simplify(A[1])
A2_3 = simplify(A[2])

#%% obtain the forward kinematics mapping [X Y Z PHI]
f_r_3D = get_f_r(T)
f_r_3D
# %%

#%% select only [X Y PHI] components
f_r = Matrix([f_r_3D[0], f_r_3D[1], f_r_3D[2]])
# select only [X Y PHI] components
f_r

#%% 
op = f_r 
wRo = Matrix([
    [-1, 0, 0],
    [0, 0, 1],
    [0, 1, 0]
            ])
wp = wRo * op
wp

#third part of the exam ------------------------------------------------
# %% declare q1, q2, q3, q4 as symbols and l1, l2 as symbols
q1, q2, q3, q4 = symbols('q1 q2 q3 q4')
l1, l2 = symbols('l1 l2')

# %% create the DH table
dh_table = Matrix([
    [0, l1, 0, q1],           # alpha, a, d, theta for joint 1
    [-pi/2, 0, 0, q2-pi/2],    # alpha, a, d, theta for joint 2
    [pi/2, 0, q3, 0],
    [0, l2, 0, q4+pi/2]
])

# %% compute transformation matrices
T, A = dh_matrix(dh_table)
T = simplify(T)
A0_1 = simplify(A[0])
A1_2 = simplify(A[1])
A2_3 = simplify(A[2])
A0_1
A1_2
A2_3

# %% get the forward kinematics mapping [X Y Z PHI]
f_r_3D = get_f_r(T)
f_r_3D

# %% select only [X Y PHI] components
f_r = Matrix([f_r_3D[0], f_r_3D[1], f_r_3D[3]])
# %% compute the jacobian
vars_list = Matrix([q1, q2, q3, q4])
j = f_r.jacobian(vars_list)
j = simplify(j)
j

# %% transform the jacobian to frame 1
R0_1 = get_rotation_mat(A[0])
j1_3 = R0_1.transpose() * j
j1_3 = simplify(j1_3)
j1_3

#convert j1_3 to sympy matrix
j1_3 = Matrix(j1_3)
#%% get minor of j1_3 using get_minors
minors = get_minors(j1_3,3)
# %% compute the determinant of the minors
determinants = [m.det() for m in minors]

#%%find the values of q1, q2, q3, q4 that make the first minor determinant zero
singular_1 = solve(determinants[1], [q1,q2,q3,q4], dict=True)
singular_1
# %%
