#%% Import required libraries
from sympy import Matrix, symbols, simplify, pi, sqrt, cos, sin, atan2
import numpy as np
from robotics_lib.geometric import euler_rotation

#%% Calculate angles
b = float(atan2(-np.sqrt(3)/2, 0.5))
a = float(atan2(0, -0.5/cos(b)))
c = float(atan2(0.5/cos(b), 0))

print(f"a = {np.degrees(a):.4f}°")
print(f"b = {np.degrees(b):.4f}°")
print(f"c = {np.degrees(c):.4f}°")

#%% Verify solution
# Create rotation matrix with computed angles
R_check = euler_rotation('yxz', [a, b, c])

# Original matrix
R = Matrix([
    [0, 1, 0],
    [0.5, 0, sqrt(3)/2],
    [sqrt(3)/2, 0, -0.5]
])

# Check difference
diff = simplify(R - R_check)
print("\nVerification - difference should be zero matrix:")
print(diff)

#%% Example singular case
R_singular = euler_rotation('yxz', [pi/2, pi/2, 0])
print("\nSingular configuration example:")
print(simplify(R_singular))