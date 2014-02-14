Australian SPF operating model and management strategy evaluation
==========================================================================================

This repository stores the most up-to-date code for the operating model used for the Australian SPF MSE. The model is available under MIT License. The model was developed by F. Hurtado-Ferro and A.E. Punt. Please reference the authors and this repository if you use this model.

Using the model
-----------------
The model is provided as an executable file (`SPFOM.exe`). To use the model, three input files are required: `SPFspp.SPEC`, `SPFOM.SPEC`, and `[xxx].DAT` (where [xxx] is the name of the species data file). `SPFspp.SPEC` tells the model which species and stock data file to use. 'SPFOM.SPEC' controls the model, harvest control rule and sensitivity scenarios. `[xxx].DAT` files store the parameters relevant to each stock. To run the model, all three input files and the executable need to be in the same folder.

The model produces an output file called `SUMMARY.OUT`. Rows are simulation years, columns are the different output quantities, organized by simulation. The first column gives the simulation year. Output quantities are organized in groups of nine columns: Biomass 1+, SSB, Recruitment, Catch, True harvest rate, Observed harvest rate, Environmental variable (not currently used), Mean age of the population, and Mean age of the catch. 

Note that only a Windows version of the model is currently available. If you need a different version, the source code for the model is provided, so you can use your favorite Fortran compiler. 

Compiling the model
---------------------
If you need to make changes to the source file, or want to compile the model to be used on a different operating system, you will need three files: `SPFOM.FOR`, `SPFOM.INC`, and `COMMON.FOR`. `SPFOM.FOR` is the main source file for the model. `SPFOM.INC` stores the variable declarations. `COMMON.FOR` stores many functions for random number generation used in the main code. 

The current version of the model was compiled using [GFortran] (http://gcc.gnu.org/wiki/GFortran), part of [GCC] (http://gcc.gnu.org/) (GNU Compiler Collection).

Other resources provided
--------------------------
`Sort.exe` reads the `SUMMARY.OUT` file and returns a new file (`PerfInd.OUT`) including mean, StdDev, and five quantiles (0.05, 0.25, 0.50, 0.75 and 0.95) for six quantities, organized by column: Biomass 1+, SSB, Catch for the whole time series, and for the last five years of the simulation. Source code for this program is also provided.

`F vs Depletion.R` is R code to produce profiles of harvest rate vs. depletion for all stocks. The code is self-contained, but you will need to specify pathways to where the model is located and where to store results. `F vs Depletion-Parallel.R` is similar, but allows for the use of parallelization to speed up runtimes. 
