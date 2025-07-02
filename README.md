# GliderSim

Simulation of equations presented in the paper, "Model-Based Feedback Control of Autonomous Underwater Gliders" by Naomi Ehrich Leonard and Joshua G. Graver (DOI: https://doi.org/10.1109/48.972106).

The program uses the SymPy, NumPy and Pandas libraries of Python for computation, substitution and storing the final results.

The program simulates the equations in the vertical (X-Z) plane for angle of attack and the position of movable mass as per the equations given on page: 639 and 640. The values of drag and lift coefficients are taken as it is from the paper along with the values of mass and mass coefficients.

A seperate program is present to calculate parameters for wings. Aspect Ratio of the wing shape is taken as input here

The inputs for the program include Glide Angle (zeta_d) and the Lift and Drag Coefficient (K_L, K_D0, K_D)






