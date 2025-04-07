# Instructions for GitHub Copilot

## Coding Standards
- Follow PEP 8 guidelines for Python code.
- Use meaningful variable names.
- Write comments for complex logic.
- Use type hints for all functions.
- Ensure all functions have docstrings explaining their purpose, parameters, and return values.

## Project-Specific Guidelines
- Use the `sympy` library for symbolic mathematics.
- Use the `numpy` library for numerical operations.
- Use `matplotlib` for plotting and visualizations.
- Use `robotics_lib` for kinematics and geometric functions.
- Use `simplify` from `sympy` to simplify symbolic expressions.
- Use `subs` from `sympy` to substitute numerical values into symbolic expressions.
- keep the same style of the examples usage in the library functions.

## Import Required Libraries, import more if needed
```python
from sympy import symbols, pi, Matrix, simplify
from robotics_lib.kinematics import dh_matrix, get_f_r, get_rotation_mat, geometric_jacobian
import numpy as np
import matplotlib.pyplot as plt
#%matplotlib inline
```