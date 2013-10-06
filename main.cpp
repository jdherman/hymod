/*
Copyright (C) 2010-2012 Josh Kollat, Jon Herman, Patrick Reed and others.

Rainfall-Runoff Models is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The Rainfall-Runoff Models is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the Rainfall-Runoff Models.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "HyMod.h"

struct HyMod hymod;

// Time period: 10/1/1961 to 9/30/1972 (1 year of warmup plus 10-year period)
int dayStartIndex = 274; //10-1 is day 274 of the year
int nDays = 4017; // length of simulation, including leap years
int startingIndex = 5023-1; //The starting index of the data file corresponding with the start date

void init_hymod(string dataFile)
{
    hymod.data.nDays = nDays;
    hymod.parameters.Nq = 3; // number of quickflow reservoirs
    hymod.parameters.Kv = 1.0; // vegetation parameter

    readMOPEXData(&hymod.data, dataFile);

    //Calculate the Hamon Potential Evaporation for the time series
    calculateHamonPE(startingIndex, nDays, dayStartIndex);

    // Allocate time series for model states
    hymod.states.snow_store     = new double [nDays];
    hymod.states.XHuz = new double [nDays];
    hymod.states.XCuz = new double [nDays];
    hymod.states.Xs   = new double [nDays];
    hymod.states.Xq   = new double* [nDays];
    for (int i=0; i < nDays; i++) 
        hymod.states.Xq[i] = new double[hymod.parameters.Nq];

    // Allocate time series for model fluxes
    hymod.fluxes.snow           = new double [nDays];
    hymod.fluxes.melt           = new double [nDays];
    hymod.fluxes.effPrecip      = new double [nDays];
    hymod.fluxes.AE             = new double [nDays];
    hymod.fluxes.OV             = new double [nDays];
    hymod.fluxes.Qq             = new double [nDays];
    hymod.fluxes.Qs             = new double [nDays];
    hymod.fluxes.Q              = new double [nDays];

    return;
}

//This is the function that gets called whenever an individual gets evaluated
void calc_hymod(double* parameters)
{
    // assign parameter values for this run
    hymod.parameters.Ks    = parameters[0]; //1.0 parameters[4]; //Ks is now specified in days
    hymod.parameters.Kq    = parameters[1]; //1.0 parameters[3]; //Kq is now specified in days
    hymod.parameters.DDF   = parameters[2];
    hymod.parameters.Tb    = parameters[3];
    hymod.parameters.Tth   = parameters[4];
    hymod.parameters.alpha = parameters[5];
    hymod.parameters.B     = parameters[6];
    hymod.parameters.Huz   = parameters[7];
    hymod.parameters.Cpar = hymod.parameters.Huz / (1.0 + hymod.parameters.B); // max capacity of soil moisture tank

    // Reinitialize everything to zero
    zero_states_and_fluxes(nDays);

    //Run Model for Simulation Period
    int dataDay;
    for (int modelDay = 0; modelDay < nDays; modelDay++)
    {
        //Since used as an index, we need to convert to zero indexing
        dataDay = startingIndex + modelDay;

        // Run snow model to find effective precip for this timestep
        hymod.fluxes.effPrecip[modelDay] = snowDD(modelDay, dataDay);

        // Run Pdm soil moisture accounting including evapotranspiration
        PDM_soil_moisture(modelDay, dataDay);

        // Run Nash Cascade routing of quickflow component
        double new_quickflow = hymod.parameters.alpha * hymod.fluxes.OV[modelDay];
        hymod.fluxes.Qq[modelDay] = Nash(hymod.parameters.Kq, hymod.parameters.Nq, new_quickflow, hymod.states.Xq[modelDay]);

        // Run Nash Cascade routing of slowflow component
        double new_slowflow = (1.0-hymod.parameters.alpha) * hymod.fluxes.OV[modelDay];
        hymod.fluxes.Qs[modelDay] = Nash(hymod.parameters.Ks, 1, new_slowflow, &hymod.states.Xs[modelDay]);

        // Set the intial states of the next time step to those of the current time step
        if (modelDay < nDays-1)
        {
            hymod.states.XHuz[modelDay+1] = hymod.states.XHuz[modelDay];
            hymod.states.Xs[modelDay+1]   = hymod.states.Xs[modelDay];
            hymod.states.snow_store[modelDay+1] = hymod.states.snow_store[modelDay];

            for(int m = 0; m < hymod.parameters.Nq; m++)
                hymod.states.Xq[modelDay+1][m] = hymod.states.Xq[modelDay][m];
        }

        hymod.fluxes.Q[modelDay] = hymod.fluxes.Qq[modelDay] + hymod.fluxes.Qs[modelDay];
    }

    return;
}

void finalize_hymod()
{
    delete[] hymod.fluxes.snow;
    delete[] hymod.fluxes.melt;
    delete[] hymod.states.snow_store;
    
    delete[] hymod.fluxes.effPrecip;
    delete[] hymod.states.XHuz;
    delete[] hymod.states.XCuz;
    delete[] hymod.states.Xs;
    delete[] hymod.fluxes.AE;
    delete[] hymod.fluxes.OV;
    delete[] hymod.fluxes.Qq;
    delete[] hymod.fluxes.Qs;
    delete[] hymod.fluxes.Q;
    
    for (int i=0; i < nDays; i++) delete[] hymod.states.Xq[i];
	delete[] hymod.states.Xq;
}

int main(int argc, char **argv)
{    
    int nParams = 8;
    double* parameters = new double [nParams]; 
    double sumQobs, sumQsim, sumPrecip;

    // initialize -- the argument is the path to the data file
    init_hymod(argv[1]);

    while(!cin.eof())
    {
        // Read model parameters from stdin
        for (int i=0; i < nParams; i++) {
			cin >> parameters[i];	
        }
		
        // Run model with this set of parameters
        calc_hymod(parameters);

        // Check observed and simulated water balance, and print results to stdout (replace with desired objectives)
        sumQobs = 0; sumQsim = 0; sumPrecip = 0;
        for (int i=0; i<nDays; i++) {
            sumQsim += hymod.fluxes.Q[i];
            sumQobs += hymod.data.flow[i];
            sumPrecip += hymod.data.precip[i];
        }

        if(!cin.eof())
            cout << "Observed: " << sumQobs << ", Simulated: " << sumQsim << ", Precip: " << sumPrecip << endl;
    }

    finalize_hymod();
    delete[] parameters;

    return 0;
}