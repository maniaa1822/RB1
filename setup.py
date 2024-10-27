from setuptools import setup, find_packages

setup(
    name="robotics_lib",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        'numpy',
        'sympy',
        'matplotlib',
        'scipy'
    ],
    author="Your Name",
    author_email="your.email@example.com",
    description="A robotics library for kinematics and dynamics",
    keywords="robotics, kinematics, dynamics",
)