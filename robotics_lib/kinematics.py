"""
robotics_lib/kinematics.py - Core robotics functions for kinematics calculations
"""

from sympy import Matrix, simplify, cos, sin, eye, zeros

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

def geometric_jacobian(f_r, sequence, q_in, dh_params):
    """
    Calculate the geometric Jacobian for a robot manipulator.
    
    Parameters:
    f_r: sympy Matrix
        Forward kinematics mapping [x, y, z, alpha_z]
    sequence: str
        A string containing the sequence of 'r's and 'p's for revolute and prismatic joints
    q_in: sympy Matrix
        Column vector of joint variables
    dh_params: sympy Matrix
        DH parameters table
        
    Returns:
    Tuple[Matrix, Matrix]:
        Jl: Linear part of the geometric Jacobian (3×n)
        Ja: Angular part of the geometric Jacobian (3×n)
    """
    # Validate inputs
    if len(sequence) != len(q_in):
        raise ValueError("Sequence length must match number of joint variables")
    
    sequence = sequence.lower()
    n_dof = len(sequence)
    
    # Calculate transformation matrices
    _, A = dh_matrix(dh_params)
    
    # Initialize Jacobians
    # Linear part (from position part of f_r)
    Jl = f_r[:3, :].jacobian(q_in)
    
    # Angular part
    Ja = zeros(3, n_dof)
    
    # Calculate cumulative rotation matrices
    R_current = eye(3)
    R_list = []  # Store intermediate rotation matrices
    
    for i in range(n_dof):
        if i == 0:
            R_list.append(get_rotation_mat(A[0]))
        else:
            R_list.append(R_list[i-1] * get_rotation_mat(A[i]))
    
    # Fill the angular part
    for i in range(n_dof):
        if sequence[i] not in ['r', 'p']:
            raise ValueError(f"Invalid joint type at position {i}. Must be 'r' or 'p'")
        
        if sequence[i] == 'r':
            # For revolute joints, use z-axis of previous frame
            if i == 0:
                Ja[:, i] = Matrix([0, 0, 1])
            else:
                prev_rot = R_list[i-1]
                Ja[:, i] = prev_rot * Matrix([0, 0, 1])
        
        # For prismatic joints, the angular part remains zero
    
    return simplify(Jl), simplify(Ja)