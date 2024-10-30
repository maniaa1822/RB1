"""
robotics_lib/__init__.py - Initialize robotics library
"""

from .kinematics import (
    dh_matrix,
    get_f_r,
    get_rotation_mat,
    geometric_jacobian
)

from .geometric import *

__all__ = [
    # Kinematics functions
    'dh_matrix',
    'get_f_r',
    'get_rotation_mat',
    'geometric_jacobian',
    
    # Geometric functions
    'elem_rot_mat',
    'elem_rot_mat_inverse',
    'euler_rotation',
    'euler_rotation_inverse',
    'rpy_rotation',
    'angle_axis_to_rot_mat',
    'angle_axis_from_rot_mat',
    'affine_T',
    'affine_get_translation',
    'affine_get_R',
    'get_minors'
]

__version__ = '0.1.0'