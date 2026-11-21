# Numerical reproduction, no observational fitting
Requirements: Python >=3.10; NumPy, SciPy, mpmath, matplotlib.
From the root of this research folder:
    python -m pip install -r code/requirements.txt
    python code/reproduce_profiles.py

The script implements the potential formulation, checks its critical-center slope, checks the published F583-4 *mass formula* from rounded printed parameters, and generates two internal-equilibrium figures.
Important distinction: the profile ODE is **not** a fit to the 22 LSB galaxy velocities, and the code does **not** provide a cosmological transfer-function result. To test the main hypothesis one must load original rotation curve measurements plus baryonic uncertainties and eventually linear cosmological structure data.
