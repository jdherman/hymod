Rainfall-Runoff Models
==========================
Hymod, HBV, and SAC-SMA Models in C/C++
--------------------------
Copyright (C) 2010-2012 Josh Kollat, Jon Herman, Patrick Reed and others. Licensed under the GNU Lesser General Public License.

Citation:
Herman, J.D., P.M. Reed, and T. Wagener (2013), Time-varying sensitivity analysis clarifies the effects of watershed model formulation on model behavior, Water Resour. Res., 49, doi:10.1002/wrcr.20124.
([Link to Paper](http://onlinelibrary.wiley.com/doi/10.1002/wrcr.20124/abstract))

Contents:
* `hymod_src/`: Hymod, an 8-parameter model based on the Probability Distributed Model (PDM) ([Moore 2007](http://hal.archives-ouvertes.fr/hal-00305633/))
* `hbv_src/`: HBV, an 11-parameter model with a similar soil moisture structure to Hymod ([Bergstrom 1995](http://www.cabdirect.org/abstracts/19961904773.html))
* `hlrms_src/`: The lumped Sacramento Soil Moisture Accounting (SAC-SMA) model, which contains 17 parameters ([Burnash and Singh 1995] (http://www.cabdirect.org/abstracts/19961904770.html))
* `shared_src/`: Files common between all three models, including objective calculations
* `example_data/`: Example precipitation and streamflow data for the Guadalupe River, Texas, from the MOPEX dataset. File paths will need to be modified accordingly in the code. Hymod and HBV use the file `GUA.in`, while SAC-SMA depends on all three files in the `GUA/` subdirectory.

In the current configuration, the model reads in a set of parameter values from a sample file. It then calculates a set of performance objectives for each solution and prints them to separate files. The objectives printed include root mean squared error (RMSE), a log-transform of RMSE (TRMSE), runoff coefficient error (ROCE), and slope of the flow duration curve error (SFDCE). These objectives are printed at the decadal, annual, seasonal, monthly, and daily timescales.

To compile and run:

* All three models can be compiled by running their respective makefiles: `cd hbv_src && make`, or `cd hymod_src && make`, etc.
* To run Hymod: `./hymod none <mopex_file> <sample_file> <num_samples> <output_file_name>`
* To run HBV: `./hbv none <mopex_file> <sample_file> <num_samples> <output_file_name>`
* To run SAC-SMA: `./hlrms <sample_file> <num_samples> <output_file_name> <mopex_file_prefix>`

Arguments:
* `<mopex_file>` in this example is the file `GUA.in`, with the appropriate path
* `<sample_file>` is a file containing parameter samples for the model
* `<num_samples>` is the number of rows in `<sample_file>`
* `<output_file_name>` is the prefix given to the set of output files that will be generated

Rainfall-Runoff Models is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Rainfall-Runoff Models is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the Rainfall-Runoff Models.  If not, see <http://www.gnu.org/licenses/>.