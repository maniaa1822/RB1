"""
robotics_lib/kinematics.py - Core robotics functions for kinematics calculations
"""

from sympy import Matrix, simplify, cos, sin, eye

def dh_matrix(dh_table):
    """
    Calculate the DH transformation matrix for robot kinematics.
    
    Parameters:
    dh_table: sympy Matrix of shape (n,4)
           Each row contains [alpha, a, d, theta] parameters
           
    Returns:
    T: The final transformation matrix (product of all individual transformations)
    A: List of individual transformation matrices
    """
    T = eye(4)
    A = []
    
    n = dh_table.shape[0]
    
    for i in range(n):
        # Extract DH parameters from the row
        line = [dh_table[i, j] for j in range(4)]
        alpha, a, d, theta = line
        
        # Create transformation matrix
        R = Matrix([
            [cos(theta), -cos(alpha)*sin(theta), sin(alpha)*sin(theta), a*cos(theta)],
            [sin(theta), cos(alpha)*cos(theta), -sin(alpha)*cos(theta), a*sin(theta)],
            [0, sin(alpha), cos(alpha), d],
            [0, 0, 0, 1]
        ])
        
        A.append(R)
        T = T * R
    
    return simplify(T), A

def get_f_r(T):
    """
    Get the forward kinematics mapping from joint space to Cartesian space.
    
    Parameters:
    T: sympy Matrix
       Homogeneous transformation matrix (symbolic)
       
    Returns:
    f_r: sympy Matrix
         Forward kinematics mapping [x, y, z, alpha_z]
    """
    # Calculate alpha_z as sum of symbolic variables in rotation matrix
    rotation_vars = set()
    for i in range(3):
        for j in range(3):
            if hasattr(T[i,j], 'free_symbols'):
                rotation_vars.update(T[i,j].free_symbols)
    
    alpha_z = sum(rotation_vars) if rotation_vars else 0
    
    # Create forward kinematics mapping
    f_r = Matrix([
        T[0,3],  # x
        T[1,3],  # y
        T[2,3],  # z
        alpha_z  # rotation
    ])
    
    return simplify(f_r)

def get_rotation_mat(T):
    """
    Extract the rotation matrix from a homogeneous transformation matrix.
    
    Parameters:
    T: sympy Matrix (4x4)
       Homogeneous transformation matrix
       
    Returns:
    R: sympy Matrix (3x3)
       Rotation matrix
    """
    return Matrix([[T[i,j] for j in range(3)] for i in range(3)])