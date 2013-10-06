###Hymod Rainfall-Runoff Model (C/C++)

###Note: This is a work in progress. It will not be ready until this note has been removed.

Hymod Rainfall-Runoff Model, based on the Probability-Distributed Model concept ([Moore 2007](http://hal.archives-ouvertes.fr/hal-00305633/)). Runs on a daily timestep and saves all states and fluxes from each day for further analysis. Currently configured to read multiple parameter sets from `stdin` and evaluate them in order, but this can be easily modified for a different application.

The model is mostly written in C (with structs instead of classes, for example), but uses a few C++ features for I/O. It has been tested for conservation of mass, and Valgrind-ed (Valground?) for memory leaks.

Contents:
* `MOPEXData.cpp/h`: Read and store forcing data from the MOPEX dataset using the format shown in the `example_data` directory. This will not be needed for users who have their own forcing data in a different format.
* `HyMod.h`: Defines the global `hymod` structure to store all states and fluxes at each timestep over the course of the evaluation.
* `HyMod.cpp`: Defines the functions for the processes in the model: degree-day snow, PDM soil moisture, Hamon PE, and the Nash cascade for the quickflow reservoirs. 
* `main.cpp`: Defines the initialization function (called once), the calculation function (called for each model evaluation), and the main function (performs model runs for each parameter set read from `stdin`).

To compile and run:

* Run `make` to compile. Modify the makefile first to use a different compiler or flags.
* Run `./hymod my_forcing_data.txt < my_parameter_samples.txt`

Arguments:
* `my_forcing_data.txt`: see the `example_data` directory for the format being used
* `my_parameter_samples.txt`: parameter sets to be evaluated in the model, with one parameter per column. Currently there are 8 parameters being read into the model, which would correspond to 8 columns per row of this file. The parameters are read from `stdin`, hence the `<` operator to pipe the contents of the file to the executable. 

In its current form, the model will output the total sum of Qobs, Qsimulated, and Precip. from the time period. However, the output can easily be modified to include any combination of states/fluxes from any time during the simulation.

Based on work from the following paper:
Herman, J.D., P.M. Reed, and T. Wagener (2013), Time-varying sensitivity analysis clarifies the effects of watershed model formulation on model behavior, Water Resour. Res., 49, doi:10.1002/wrcr.20124.
([Link to Paper](http://onlinelibrary.wiley.com/doi/10.1002/wrcr.20124/abstract))

Copyright (C) 2010-2013 Jon Herman, Josh Kollat, and others.

Hymod is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Hymod is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with Hymod.  If not, see <http://www.gnu.org/licenses/>.
