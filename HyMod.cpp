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

void zero_states_and_fluxes(int ndays)
{
    for (int k=0; k < ndays; k++)
    {
        hymod.fluxes.snow[k]  = 0.0;
        hymod.fluxes.melt[k]  = 0.0;
        hymod.fluxes.effPrecip[k] = 0.0;
        hymod.fluxes.AE[k] = 0.0;
        hymod.fluxes.OV[k] = 0.0;
        hymod.fluxes.Qq[k] = 0.0;
        hymod.fluxes.Qs[k] = 0.0;
        hymod.fluxes.Q[k]  = 0.0;

        hymod.states.snow_store[k] = 0.0;
        hymod.states.XHuz[k]      = 0.0;
        hymod.states.XCuz[k] = 0.0;
        hymod.states.Xs[k] = 0.0;
        for (int m=0; m < hymod.parameters.Nq; m++) hymod.states.Xq[k][m] = 0.0;
    }
    return;
}

// Allocate time series arrays for states, fluxes, and forcing data
void hymod_allocate(int ndays)
{
    hymod.states.snow_store     = new double [ndays];
    hymod.states.XHuz = new double [ndays];
    hymod.states.XCuz = new double [ndays];
    hymod.states.Xs   = new double [ndays];
    hymod.states.Xq   = new double* [ndays];
    for (int i=0; i < ndays; i++) 
        hymod.states.Xq[i] = new double[hymod.parameters.Nq];

    hymod.fluxes.snow           = new double [ndays];
    hymod.fluxes.melt           = new double [ndays];
    hymod.fluxes.effPrecip      = new double [ndays];
    hymod.fluxes.AE             = new double [ndays];
    hymod.fluxes.OV             = new double [ndays];
    hymod.fluxes.Qq             = new double [ndays];
    hymod.fluxes.Qs             = new double [ndays];
    hymod.fluxes.Q              = new double [ndays];
    hymod.fluxes.PE = new double [ndays];
}

// clean up memory
void hymod_delete(int ndays) 
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
    
    for (int i=0; i < ndays; i++) delete[] hymod.states.Xq[i];
    delete[] hymod.states.Xq;
}

void PDM_soil_moisture(int modelDay, int dataDay)
{
    double Cbeg, OV2, PPinf, Hint, Cint, OV1; // temporary variables for intermediate calculations
    
    // Storage contents at begining
    Cbeg = hymod.parameters.Cpar * (1.0 - pow(1.0-(hymod.states.XHuz[modelDay]/hymod.parameters.Huz),1.0+hymod.parameters.B));

    // Compute overflow from soil moisture storage element
    OV2 = max(0.0, hymod.fluxes.effPrecip[modelDay] + hymod.states.XHuz[modelDay] - hymod.parameters.Huz);

    // Remaining net rainfall
    PPinf = hymod.fluxes.effPrecip[modelDay] - OV2;

    // New actual height in the soil moisture storage element
    Hint = min(hymod.parameters.Huz, hymod.states.XHuz[modelDay] + PPinf);

    // New storage content
    Cint = hymod.parameters.Cpar*(1.0-pow(1.0-(Hint/hymod.parameters.Huz),1.0+hymod.parameters.B));

    // Additional effective rainfall produced by overflow from stores smaller than Cmax
    OV1 = max(0.0, PPinf + Cbeg - Cint);

    // Compute total overflow from soil moisture storage element
    hymod.fluxes.OV[modelDay] = OV1 + OV2;
    
    // Compute actual evapotranspiration
    hymod.fluxes.AE[modelDay] = min(Cint, (Cint/hymod.parameters.Cpar)*hymod.fluxes.PE[modelDay]*hymod.parameters.Kv);
    
    // Storage contents and height after ET occurs
    hymod.states.XCuz[modelDay] = max(0.0, Cint - hymod.fluxes.AE[modelDay]);
    hymod.states.XHuz[modelDay] = hymod.parameters.Huz*(1.0-pow(1.0-(hymod.states.XCuz[modelDay]/hymod.parameters.Cpar),1.0/(1.0+hymod.parameters.B)));

    return;

}

double Nash(double K, int N, double Qin, double *X)
{
    //Initialization
    double *OO = new double[N];
    double Qout;                       //Flow out of series of reservoirs
    
    //Loop through reservoirs
    for (int Res = 0; Res < N; Res++)
    {
        OO[Res] = K*X[Res];
        X[Res]  = X[Res] - OO[Res];

        if (Res==0) X[Res] = X[Res] + Qin; 
        else        X[Res] = X[Res] + OO[Res-1];
    }

    // The outflow from the cascade is the outflow from the last reservoir
    Qout = OO[N-1];
    delete[] OO;
    return Qout;
}

double snowDD(int modelDay, int dataDay)
{
    double Qout; // effective precip after freezing/melting

    //If temperature is lower than threshold, precip is all snow
    if (hymod.data.avgTemp[dataDay] < hymod.parameters.Tth)
    {
        hymod.fluxes.snow[modelDay] = hymod.data.precip[dataDay];
        Qout = 0.0;
    }
    else //Otherwise, there is no snow and it's all rain
    {
        hymod.fluxes.snow[modelDay] = 0.0;
        Qout = hymod.data.precip[dataDay];
    }

    //Add to the snow storage for this day
    hymod.states.snow_store[modelDay] += hymod.fluxes.snow[modelDay];

    //Snow melt occurs if we are above the base temperature (either a fraction of the store, or the whole thing)
    if (hymod.data.avgTemp[dataDay] > hymod.parameters.Tb)
    {
        hymod.fluxes.melt[modelDay] = min(hymod.parameters.DDF*(hymod.data.avgTemp[dataDay]-hymod.parameters.Tb), hymod.states.snow_store[modelDay]);
    }
    //Otherwise, snowmelt is zero
    else
    {
        hymod.fluxes.melt[modelDay] = 0.0;
    }

    //Update the snow storage depending on melt
    hymod.states.snow_store[modelDay] -= hymod.fluxes.melt[modelDay];
    if(hymod.states.snow_store[modelDay] < 0.0) hymod.states.snow_store[modelDay] = 0.0;

    //Qout is any rain + snow melt
    Qout += hymod.fluxes.melt[modelDay];

    return Qout;
}

void calculateHamonPE(int dataIndex, int nDays, int startDay)
{
    int oldYear;
    int counter;
    double evap_P, evap_day_length, evap_eStar;

    //Initialize the starting year
    oldYear = hymod.data.date[dataIndex][0];
    counter = startDay-1;

    //Fill out each of the arrays
    for (int i=0; i<nDays; i++)
    {
        //If the years hasn't changed, increment counter
        if (hymod.data.date[dataIndex+i][0] == oldYear) counter++;
        //If it has changed, reset counter - this handles leap years
        else counter = 1;

        evap_P = asin(0.39795*cos(0.2163108 + 2.0 * atan(0.9671396*tan(0.00860*double(counter-186)))));
        evap_day_length = 24.0 - (24.0/PI)*(acos((sin(0.8333*PI/180.0)+sin(hymod.data.gageLat*PI/180.0)*sin(evap_P))/(cos(hymod.data.gageLat*PI/180.0)*cos(evap_P))));
        evap_eStar = 0.6108*exp((17.27*hymod.data.avgTemp[dataIndex+i])/(237.3+hymod.data.avgTemp[dataIndex+i]));
        hymod.fluxes.PE[i] = (715.5*evap_day_length*evap_eStar/24.0)/(hymod.data.avgTemp[dataIndex+i] + 273.2);

        oldYear = hymod.data.date[dataIndex+i][0];
    }

    return;
}