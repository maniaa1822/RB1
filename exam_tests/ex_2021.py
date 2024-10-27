# rotation with euler angles ex 2021

#%%import 
from robotics_lib.geometric import *
from robotics_lib.kinematics import *
from sympy import symbols, Matrix, pi, cos, sin, simplify, eye, sqrt,deg

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

