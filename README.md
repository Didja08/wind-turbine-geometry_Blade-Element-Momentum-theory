# wind-turbine-geometry_Blade-Element-Momentum-theory
This MATLAB code defines a wind turbine geometry and calculates various aerodynamic parameters using Blade Element Momentum theory
The first section of the code defines the geometry of the wind turbine blades using the rotor radius, number of blades, number of blade sections, and radial and angular coordinates. The flow conditions are then defined, including the wind speed, air density, and air viscosity.

The code then interpolates airfoil properties for each blade section using the angle of attack, lift coefficient, and drag coefficient. Lift and drag forces are calculated for each blade section using the interpolated airfoil properties, and the normal and tangential force coefficients are calculated. Axial and tangential induction factors are calculated, and the incremental inflow angle is determined. The pressure, viscous, and momentum forces are then calculated. Finally, the sound pressure field is computed.

This code is useful for simulating the aerodynamic performance of a wind turbine blade and can help in the design process.

##################################################################################################################################################################

To run this code, you will need to provide the following data:

N - number of blade sections
R - radius of the turbine
Omega - rotational speed of the turbine
V - wind speed
chord - chord length of each blade section
twist - twist angle of each blade section
airfoil - airfoil data (lift and drag coefficients vs. angle of attack)
rho - density of air
mu - dynamic viscosity of air
The data can be provided as arrays or matrices in the MATLAB workspace before running the code.
