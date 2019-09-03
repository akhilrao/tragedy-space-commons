# tragedy-space-commons

This repository contains data, code, and images for the paper "Orbital-use fees could more than quadruple the value of the space industry". The output was generated in R version 3.6.1.

File structure:
	- /bin contains all of the R code
	- /data contains all of the calibration and output data
	- /images contains all generated figures

To run the code and generate the paper's figures and data:
	1. Make sure R and the necessary R packages are installed. 
	You can find a list of necessary packages in the header of "/data/final_script.r"
	2. Run "final_script.r" (_R code: source("final_script.r")_)
	Note that the default is to use 3 cores. You will want to adjust this in "final_script.r".

