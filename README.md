# monte_carlo_simulation
MATLAB code for modelling of a agglomeration process in spray fluidized bed.

# Description of the code files

Monte_Carlo_milkpowder.m - main code file runs the main Monte Carlo simulation for milk powder agglomeration in spray fluidized bed.
Function_diameter.m- This function calculates the average particle size of the agglomerates, such as d32 and d50 using report A.
Function_Frequency.m - this function calculates the collision frequency for the time interval in the main code.
Function_height.m - this function calculates the imbibition height of the droplet deposition on milk powder particle surface.
Analyzing_all_growthrate.m - this code is used to analyze the growth of the agglomerate with time for all the simulation that runs for different operating condition

# How to run the code

To start the simulation just run the monte_Carlo_milkpowder.m file in Matlab.

During the simulation the relative mean agglomerate diameter are shown in the command window and simulation will finish when the time screens reach at desired value, once it's finished the results can be analyzed by giving the complete name and date, for example, .'CVMC_Real_SMP_final-11-Jul-2023-reportA-' to the Analyzing_all_growthrate.m program
