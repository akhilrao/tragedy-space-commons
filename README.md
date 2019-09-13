# tragedy-space-commons

This repository contains data, code, and images for the paper "Orbital-use fees could more than quadruple the value of the space industry". The majority of the model output was generated in R 3.6.1, but "/data/bootstrap_sims.csv" with 250 bootstrapped paths was generated in R 3.4.4. Changes to the behavior of RNG in 3.6.0 mean that replicating the 250 draws will require adding "RNGkind(sample.kind="Rounding")" to the top of "/data/final_script.r".

File structure:
	- /bin contains all of the R code
	- /data contains all of the calibration and output data
	- /images contains all generated figures

To run the code and generate the paper's figures and data:
	1. Make sure R and the necessary R packages are installed. 
	You can find a list of necessary packages in the header of "/data/final_script.r"
	2. Run "final_script.r" (_R code: source("final_script.r")_)
	Note that the default is to use 3 cores. You will want to adjust this in "final_script.r".

