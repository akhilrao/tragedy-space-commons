# tragedy-space-commons

This repository contains data, code, and images for the paper "Orbital-use fees could more than quadruple the value of the space industry". The output was generated in R version 3.6.1.

### How is this repository organized?

	* /bin contains all of the R code
	* /data contains all of the calibration and output data
	* /images contains all generated figures

### How do I recreate the figures and data files?

To run the code and generate the paper's figures and data:

	1. Make sure R and the necessary R packages are installed. 
	You can find a list of necessary packages in the header of "/data/final_script.r"
	2. Run "final_script.r" (_R code: source("final_script.r")_)
	Note that the default is to use 3 cores. You will want to adjust this in "final_script.r".

### What does the code in /bin do?

Broadly, there are two types of code files: those which describe functions or algorithms to compute functions ("methods"), and those which contain scripts that compute those functions/implement those algorithms ("content").

There are four main "methods" files:

	* `equations.r' contains the main theoretical equations. If there's a (simple!) economic or physical equation defined in the SI, it's in here. 
	* `simulation_algorithms.r' contains the computational algorithms used and described in the SI. Dynamic programming, path generation -- it's in here, though it's broken out into many smaller functions.
	* `simulation_functions.r' contains the functions used in `simulation_algorithms.r'. These include grid generators and a linear interpolator. There's some overlap here with the two files above.
	* `plotting_functions.r' contains some functions used to generate plots, mostly in the calibration section.

The remaining files are "content" files: