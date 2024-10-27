from typing import List, Union, Tuple
import sympy as sp
from sympy import Matrix, cos, sin, simplify, eye, atan2, sqrt, acos, asin
from sympy import Matrix, zeros, ones, eye, Expr, Symbol, simplify

def elem_rot_mat(axis: str, angle: sp.Symbol) -> sp.Matrix:
    """
    Create an elementary rotation matrix about a specified axis.
    
    Parameters:
    -----------
    axis: str
        The axis about which to rotate ('x', 'y', or 'z')
    angle: sp.Symbol
        The angle of rotation (symbolic)
        
    Returns:
    --------
    R: sp.Matrix
        The elementary rotation matrix
    """
    axis = axis.lower()
    if axis == 'x':
        R = Matrix([
            [1, 0, 0],
            [0, cos(angle), -sin(angle)],
            [0, sin(angle), cos(angle)]
        ])
    elif axis == 'y':
        R = Matrix([
            [cos(angle), 0, sin(angle)],
            [0, 1, 0],
            [-sin(angle), 0, cos(angle)]
        ])
    elif axis == 'z':
        R = Matrix([
            [cos(angle), -sin(angle), 0],
            [sin(angle), cos(angle), 0],
            [0, 0, 1]
        ])
    else:
        raise ValueError("Axis must be 'x', 'y', or 'z'")
    
    return R

def euler_rotation(sequence: str, angles: List[sp.Symbol]) -> sp.Matrix:
    """
    Compute Euler rotation matrix for given sequence and angles.
    
    Parameters:
    -----------
    sequence: str
        String specifying rotation sequence (e.g., 'xyz', 'zyx')
    angles: list
        List of three symbolic angles
        
    Returns:
    --------
    R: sp.Matrix
        The resulting rotation matrix
    """
    if len(sequence) != 3 or len(angles) != 3:
        raise ValueError("Sequence and angles must both be of length 3")
    
    sequence = sequence.lower()
    
    # Check for invalid sequences
    if sequence[0] == sequence[1] or sequence[1] == sequence[2]:
        raise ValueError("Invalid Euler sequence: consecutive rotations about same axis")
    
    R = elem_rot_mat(sequence[0], angles[0])
    R = R @ elem_rot_mat(sequence[1], angles[1])
    R = R @ elem_rot_mat(sequence[2], angles[2])
    
    return simplify(R)

def rpy_rotation(sequence: str, angles: List[sp.Symbol]) -> sp.Matrix:
    """
    Compute Roll-Pitch-Yaw rotation matrix for given sequence and angles.
    RPY rotations work about fixed axes (unlike Euler rotations which work about moving axes).
    
    Parameters:
    -----------
    sequence: str
        String specifying rotation sequence (e.g., 'xyz', 'zyx')
    angles: list
        List of three symbolic angles
        
    Returns:
    --------
    R: sp.Matrix
        The resulting rotation matrix
    """
    # RPY is equivalent to reversed Euler sequence with reversed angles
    return euler_rotation(sequence[::-1], angles[::-1])

def angle_axis_rotation_direct(r: sp.Matrix, theta: sp.Symbol) -> sp.Matrix:
    """
    Compute rotation matrix from axis-angle representation (Rodrigues' formula).
    
    Parameters:
    -----------
    r: sp.Matrix
        Unit vector specifying rotation axis
    theta: sp.Symbol
        Rotation angle
        
    Returns:
    --------
    R: sp.Matrix
        The resulting rotation matrix
    """
    I = eye(3)
    S = Matrix([
        [0, -r[2], r[1]],
        [r[2], 0, -r[0]],
        [-r[1], r[0], 0]
    ])
    R = r * r.transpose() + (I - r * r.transpose()) * cos(theta) + S * sin(theta)
    return simplify(R)

def affine_T(rotation_matrix: sp.Matrix, translation: sp.Matrix) -> sp.Matrix:
    """
    Create homogeneous transformation matrix from rotation matrix and translation vector.
    
    Parameters:
    -----------
    rotation_matrix: sp.Matrix
        3x3 rotation matrix
    translation: sp.Matrix
        3x1 translation vector
        
    Returns:
    --------
    T: sp.Matrix
        4x4 homogeneous transformation matrix
    """
    T = Matrix.zeros(4, 4)
    T[:3, :3] = rotation_matrix
    T[:3, 3] = translation
    T[3, 3] = 1
    return simplify(T)

def affine_get_translation(T: sp.Matrix) -> sp.Matrix:
    """Extract translation vector from homogeneous transformation matrix."""
    return T[:3, 3]

def affine_get_R(T: sp.Matrix) -> sp.Matrix:
    """Extract rotation matrix from homogeneous transformation matrix."""
    return T[:3, :3]

def elem_rot_mat_inverse(axis: str, rotation_matrix: sp.Matrix) -> sp.Expr:
    """
    Extract the rotation angle about a given axis from a rotation matrix.
    
    Parameters:
    -----------
    axis: str
        The axis of rotation ('x', 'y', or 'z')
    rotation_matrix: sp.Matrix
        3x3 rotation matrix
        
    Returns:
    --------
    angle: sp.Expr
        The symbolic rotation angle
    """
    axis = axis.lower()
    
    if axis == 'x':
        # R = [1     0        0    ]
        #     [0   cos(θ)  -sin(θ) ]
        #     [0   sin(θ)   cos(θ) ]
        return atan2(rotation_matrix[2, 1], rotation_matrix[2, 2])
    
    elif axis == 'y':
        # R = [ cos(θ)   0   sin(θ) ]
        #     [   0      1     0    ]
        #     [-sin(θ)   0   cos(θ) ]
        return atan2(rotation_matrix[0, 2], rotation_matrix[0, 0])
    
    elif axis == 'z':
        # R = [cos(θ)  -sin(θ)   0]
        #     [sin(θ)   cos(θ)   0]
        #     [  0        0      1]
        return atan2(rotation_matrix[1, 0], rotation_matrix[0, 0])
    
    else:
        raise ValueError("Axis must be 'x', 'y', or 'z'")

def euler_rotation_inverse(sequence: str, rotation_matrix: sp.Matrix, positive_solution: bool = True) -> Tuple[sp.Expr, sp.Expr, sp.Expr]:
    """
    Extract Euler angles from a rotation matrix.
    
    Parameters:
    -----------
    sequence: str
        String specifying rotation sequence (e.g., 'xyz', 'zyx')
    rotation_matrix: sp.Matrix
        3x3 rotation matrix
    positive_solution: bool
        If True, returns the solution with positive theta when applicable
        
    Returns:
    --------
    alpha, beta, gamma: tuple of sp.Expr
        The three Euler angles
    """
    sequence = sequence.lower()
    R = rotation_matrix
    
    if len(sequence) != 3:
        raise ValueError("Sequence must be of length 3")
        
    if sequence == 'zyx':
        # Most common case for RPY
        beta = asin(-R[0,2]) if positive_solution else pi - asin(-R[0,2])
        alpha = atan2(R[1,2]/cos(beta), R[2,2]/cos(beta))
        gamma = atan2(R[0,1]/cos(beta), R[0,0]/cos(beta))
        return alpha, beta, gamma
        
    elif sequence == 'xyz':
        # Another common case
        beta = asin(R[2,0]) if positive_solution else pi - asin(R[2,0])
        alpha = atan2(-R[2,1]/cos(beta), R[2,2]/cos(beta))
        gamma = atan2(-R[1,0]/cos(beta), R[0,0]/cos(beta))
        return alpha, beta, gamma
        
    elif sequence == 'zyz':
        # Example of proper Euler angles
        beta = acos(R[2,2]) if positive_solution else -acos(R[2,2])
        alpha = atan2(R[1,2], R[0,2])
        gamma = atan2(R[2,1], -R[2,0])
        return alpha, beta, gamma
        
    else:
        raise ValueError(f"Sequence {sequence} not implemented yet")

def angle_axis_to_rot_mat(axis: sp.Matrix, angle: sp.Symbol) -> sp.Matrix:
    """
    Convert angle-axis representation to rotation matrix using Rodrigues' formula.
    
    Parameters:
    -----------
    axis: sp.Matrix
        3x1 unit vector representing rotation axis
    angle: sp.Symbol
        Rotation angle
        
    Returns:
    --------
    R: sp.Matrix
        3x3 rotation matrix
    """
    # Ensure axis is normalized
    axis = axis / sqrt(axis[0]**2 + axis[1]**2 + axis[2]**2)
    
    # Rodrigues' formula components
    K = Matrix([
        [0, -axis[2], axis[1]],
        [axis[2], 0, -axis[0]],
        [-axis[1], axis[0], 0]
    ])
    
    I = eye(3)
    R = I + sin(angle)*K + (1 - cos(angle))*K*K
    return simplify(R)

def angle_axis_from_rot_mat(R: sp.Matrix) -> Tuple[sp.Matrix, sp.Expr]:
    """
    Extract axis and angle from rotation matrix.
    
    Parameters:
    -----------
    R: sp.Matrix
        3x3 rotation matrix
        
    Returns:
    --------
    axis: sp.Matrix
        3x1 unit vector representing rotation axis
    angle: sp.Expr
        Rotation angle
    """
    # The angle can be found from the trace
    angle = acos((R.trace() - 1)/2)
    
    # The axis is related to the skew-symmetric part
    axis = Matrix([
        R[2,1] - R[1,2],
        R[0,2] - R[2,0],
        R[1,0] - R[0,1]
    ])
    
    # Normalize the axis
    axis = axis / (2*sin(angle))
    
    return simplify(axis), simplify(angle)

def get_minors(matrix: Matrix, size_minor_matrix: int) -> List[Matrix]:
    from itertools import combinations
    """
    Get all possible minors of a given size from a matrix using SymPy's built-in functions.
    
    Args:
        matrix: SymPy Matrix to extract minors from
        size_minor_matrix: Size of the minor matrices to extract (must be less than or equal to matrix dimensions)
        
    Returns:
        List of all possible minor matrices of the specified size
        
    Example:
        >>> from sympy import Matrix, symbols
        >>> x, y = symbols('x y')
        >>> M = Matrix([[x, y, 1],
        ...            [2, 3, 4],
        ...            [5, 6, 7]])
        >>> minors = get_minors(M, 2)
        >>> print(f"Number of {2}x{2} minors: {len(minors)}")
        >>> print("First minor:")
        >>> print(minors[0])
    """
    m, n = matrix.shape
    
    # Validate input
    if size_minor_matrix > min(m, n):
        raise ValueError(f"Minor size {size_minor_matrix} must be less than or equal to minimum matrix dimension {min(m, n)}")
    if size_minor_matrix <= 0:
        raise ValueError("Minor size must be positive")
        
    # Get all possible row and column combinations
    row_combinations = list(combinations(range(m), size_minor_matrix))
    col_combinations = list(combinations(range(n), size_minor_matrix))
    
    # Generate all possible minors using matrix.extract()
    minors = []
    for rows in row_combinations:
        for cols in col_combinations:
            minor = matrix.extract(rows, cols)
            minors.append(minor)
            
    return minors