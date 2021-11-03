# 4DFUCCI

4D FUCCI README File

This README file supplements the code made publicly available on GitHub. It contains a brief guide on
how the available code is designed to be run for the reader’s testing, as well as recommendations for some
parameter choices depending on available computational power.

System Requirements
===================
* MATLAB R2018a or later
* MEX and C compiler
* Additional note: High performance computing is recommended, with 4 CPUs and 32 GB RAM if setting `I = 201`. Otherwise, a personal computer is sufficient.


Installation
============
MATLAB can interpret C code with the mex function. Before running the IBM, ensure that the cells.c file is compiled with MATLAB, using the command 
```matlab 
mex cells.c -R2018a
```
in the command window in MATLAB. A Windows compatible compilation file is provided, (`cells.mexw64`), but it is recommended that the `cells.c` file is manually compiled, and required if running on Linux or MacOS, as the interpreting file is not compatible across operating systems 

Visualisation of the 3D nutrient profiles require Python, and the packages Plotly, Numpy, Scipy, and Matplotlib. Python can be installed from [here](https://www.python.org/downloads/), and the packages can be installed with the command,
```
pip install [package]
```
in the Command Prompt or Unix terminal.

Running the IBM
===============
Two approaches to running the model are provided in the repository. We make the full IBM simulation code available, 
as well as an example production of figures with a quick iteration of the simulation.

IBM Code
--------
The IBM itself is contained in the `ibm3d.m` and `cells.c` functions. These are called with parameters supplied by `hpc_script.m` and `quick_demo.m`.

Full IBM Modelling
------------------
The full IBM is recommended for use in high performance computing. To run the full suite of simulations, call the command
```matlab 
for i = 1:10
   hpc_script(i)
end
``` 
The input to `hpc_script.m`, `sim_id`, is a parameter that identifies the simulation, and sets the initial radius and population. The variability case is provided, and can be uncommented (variables `exp_N0` and `exp_rad`) to select the variable conditions. Also, `sim_id` acts as a seed for the random number generation, and names the workspaces when saving. The function `hpc_script.m` saves these workspaces into the Workspaces folder, ready for analysis and figure generation. 

Quick example
-------------
A quick example of the IBM with a reduced number of nodes (`I = 51`) and experimental period (`T = 144 h`) is also provided. This generates all figures created with the full IBM (with the exception of days 7-10 and the Python figures) and is a good, quick example realisation of the IBM. The runtime is approximately 5-6 minutes, including data analysis and figure generation. 

Reproducing Figures
===================
With 10 workspaces, the data for figures can be prepared, and the figures can be generated. All code required to reproduce the figures is readily available, however the data plotted in the figures must first be taken from the IBM and prepared for visualisation. To achieve this preparation, run the `workspace_analysis.m` file. This can be a computationally expensive process and take a large amount of time to complete. To get a quicker example of data visualisation for all 10 days of the simulation, set `runcount = 1` on line 52. The file `cMap.mat` is the workspace variable that defines the colourmap for nutrient profiles. This is loaded in automatically in the code. 

When `workspace_analysis.m` is complete, each figure can be generated by running the corresponding .m file in their corresponding folders. All figures will be saved with the day they correspond to (if appropriate) in their appropriate folders, generated by the figure creation file. For generating 3D nutrient profile plots, one 3D nutrient profile is generated by running `nutrient3Dplot.py`. To change which time point is plotted, change the variable `c_X` on line 46 of `nutrient3Dplot.py`, so that X is between 0 and 10, inclusive. This X corresponds to which day of the simulation is plotted. One nutrient profile per run of `nutrient3Dplot.py` is plotted so that the Plotly API is not overcrowded, and won't crash the tab in the browser.

Experimental Data
=================
The experimental radius data are collected in the `RadiusData` folder. `IncucyteData.csv` contains the data for external radii for a number of experiments, including the radius data for WM793b, with an initial seeding density of 10,000 cells. Files with the name format `dayX.csv` are the outer, arrested, and necrotic radius estimates from the image processing for experimental images. The data in this folder are used for the generation of Figure 8. 
