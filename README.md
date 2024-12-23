# The One-Dimensional Lateral Infiltration of Water through Porous Snow

This folder contains code to simulate the fluidic and thermodynamic evolution of water laterally flowing through a snowpack. The key variables are water depth, ice porosity, water temperature, and ice temperature.

The wrapper code is `WWCS_wrapper` which calls the function `WWCS_func_rev6`, as this project is still in progress. The function takes in water source temperature, ice sink temperature, and the desired run time to return depth, porosity, and temperature matrices, as well as some other variables to facilitate the presentation of the simulation results. The wrapper plots the results, with darker curves denoting later snapshots into the simulation. Please find some finer details about the code below.

The four key variables are returned from the function into the wrapper as m x n sized matrices, where m is the number of time steps and n is the number of spatial steps.

Physical parameters used in the simulations are described and listed with comments between lines 17-36 of the function. Certain parameters of interest include the horizontal lengthscale X and the depthscale H, which may be changed as desired.

In our nondimensionalization of the governing eqautions, we uncover a timescale T and the nondimensional groups Da (Darcy number), Ht (heat transfer timescales ratio), Pe (Péclet number), and St (Stefan number). These are used in the function and are returned to the wrapper.

To solve our governing equations, we initialize our key variables and give them realistic initial and boundary conditions. Please see the "Initial, Boundary Conditions" section of the function for the setup. We then update our mesh in the "Fill Mesh" section of the function.

The wrapper feeds the scenarios into the function, calling different source and sink temperatures and giving desired run times. For each scenario, the wrapper plots and saves the data for each of the key variables as .png and .mat files. The .png files are images of the specifically produced figures in the code while the .mat files save the entire matrices for each key variable produced from solving the governing equations.
