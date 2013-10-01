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

/*
%%=========================================================================
%% INPUTS
%%   data         = structure holding input data such as precip
%%   period       = simulation period array holding the indices into the data array
%%   periodLength = length of simulation period
%%   pars         = HyMod parameters
%%   snowPars     = Degree Day snow parameters
%%   inState      = initial states of HyMod stores
% OUTPUTS
%%   model        = model computed variables - final values
%%   snowModel    = snow model computed variables - final values
%%=========================================================================
*/

void HyModMo(MOPEXData data, HamonEvap evap, int *period, int periodLength, 
             HyModPars pars, SnowDDPars snowPars, HyModInStates inState, 
             HyModModel *model, SnowDDModel *snowModel)
{

    //Initial state values for quickflow routing tanks
    for (int i=0; i<pars.Nq; i++) model->Xq[0][i] = 0.0;

    //Initial state value for slowflow routing tank
    model->Xs[0] = 0.0;

    //Initial height corresponding to SMA tank contents
    model->XHuz[0] = 0.0;

    //Calculate b
    //pars.b = log(1-(pars.B/2.0))/log(0.5);  //We are taking out this scaling that assumes the upper limit of b is 2
    pars.b = pars.B;

    //Calculate maximum capacity of soil moisture accounting tank
    pars.Cpar = pars.Huz / (1.0+pars.b);

    //Run Model Simulation TIME Loop
    int Nday = periodLength;
    int dataDay;

    for (int modelDay=0; modelDay<Nday; modelDay++)
    {
        //Since used as an index, we need to convert to zero indexing
        //Note: the values contained in "period" are already zero based indices into the data array (i.e., don't subtract 1)
        dataDay = period[modelDay];

        //Run snow model if needed
        if (snowPars.useSnowDD == true) 
        {
            model->effPrecip[modelDay] = snowDD(modelDay, dataDay, data, snowPars, snowModel);
        }
        //If not, the effective precip is just the actual precip
        else 
        {
			model->effPrecip[modelDay] = (data.precip[dataDay]);
        }

        //Run Pdm soil moisture accounting including evapotranspiration
        Pdm03(modelDay, dataDay, pars, data, evap, model);

        //Run Nash Cascade routing of quickflow component
        model->Qq[modelDay] = Nash(pars.Kq, pars.Nq, pars.alpha*model->OV[modelDay], model->Xq[modelDay]);

        //Run Nash Cascade routing of slowflow component
        model->Qs[modelDay] = Nash(pars.Ks, 1, (1.0-pars.alpha)*model->OV[modelDay], &model->Xs[modelDay]);

        //Set the intial states of the next time step to those of the current time step
        if (modelDay<Nday-1)
        {
            //Carry the Hymod states to the next time step
            model->XHuz[modelDay+1] = model->XHuz[modelDay];
            model->Xq[modelDay+1]   = model->Xq[modelDay];
            model->Xs[modelDay+1]   = model->Xs[modelDay];
            //And be sure to carry the snow storage to the next time step
            snowModel->store[modelDay+1] = snowModel->store[modelDay];
        }
    }

    //Finalize variables   
    for (int modelDay=0; modelDay<Nday; modelDay++)
    {
        //Model computed total streamflow flux
        model->Q[modelDay] = model->Qq[modelDay] + model->Qs[modelDay];
    }

    return;

}

void Pdm03(int modelDay, int dataDay, HyModPars pars, MOPEXData data, HamonEvap evap, HyModModel *model)
{
    double Cbeg, OV2, PPinf, Hint, Cint, OV1;
    //double Sbeg, u1k, u2k, r_star, Ck, Sk;
    
    //Storage contents at begining
    Cbeg = pars.Cpar * (1.0 - pow(1.0-(model->XHuz[modelDay]/pars.Huz),1.0+pars.b));

    //Compute effective rainfall filling all storage elements
    OV2 = max(0.0, model->effPrecip[modelDay] + model->XHuz[modelDay] - pars.Huz);

    //Remaining net rainfall
    PPinf = model->effPrecip[modelDay] - OV2;

    //New actual height
    Hint = min(pars.Huz, model->XHuz[modelDay]+PPinf);

    //New storage content
    Cint = pars.Cpar*(1.0-pow(1.0-(Hint/pars.Huz),1.0+pars.b));

    //Additional effective rainfall produced by stores smaller than Cmax
    OV1 = max(0.0, PPinf + Cbeg - Cint);

    //Compute total effective rainfall
    model->OV[modelDay] = OV1 + OV2;
    
    //Computer actual evapotranspiration
    model->AE[modelDay] = min(Cint,(Cint/pars.Cpar)*evap.PE[modelDay]*pars.Kv);
    
    //Storage contents after ET
    model->XCuz[modelDay] = max(0.0,Cint - model->AE[modelDay]);
    
    //Compute final height of the reservoir
    model->XHuz[modelDay] = pars.Huz*(1.0-pow(1.0-(model->XCuz[modelDay]/pars.Cpar),1.0/(1.0+pars.b)));

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

    //Get Qout
    Qout = OO[N-1];

    //Clean up
    delete[] OO;

    //Return the flow output
    return Qout;

}

double snowDD(int modelDay, int dataDay, MOPEXData data, SnowDDPars pars, SnowDDModel *snowModel)
{

    double Qout;

    //If temperature is lower than threshold...
    if (data.avgTemp[dataDay] < pars.Tth)
    {
        //Precip is all snow
        snowModel->snow[modelDay] = data.precip[dataDay];
        //And there is no rain
        Qout = 0.0;
    }
    //Otherwise...
    else
    {
        //There is no snow
        snowModel->snow[modelDay] = 0.0;
        //And it is all rain
        Qout = data.precip[dataDay];
    }

    //Add to the snow storage for this day
    snowModel->store[modelDay] += snowModel->snow[modelDay];

    //Snow melt occurs if we are above the base temperature
    if (data.avgTemp[dataDay] > pars.Tb)
    {
        //Snow melt is either a fraction of the store, or the whole thing
        snowModel->melt[modelDay] = min(pars.DDF*(data.avgTemp[dataDay]-pars.Tb),snowModel->store[modelDay]);
    }
    //Otherwise
    else
    {
        //Snow melt is zero
        snowModel->melt[modelDay] = 0.0;
    }

    //Update the snow storage depending on melt
    snowModel->store[modelDay] -= snowModel->melt[modelDay];
    if (snowModel->store[modelDay] < 0.0) snowModel->store[modelDay] = 0.0;

    //Qout is any rain + snow melt
    Qout += snowModel->melt[modelDay];

    return Qout;
}

void calculateHamonPE(MOPEXData *data, int dataIndex, int nDays, HamonEvap *evap, int startDay)
{
    int oldYear;
    int counter;

    evap->PE = new double [nDays];

    //Initialize the starting year
    oldYear = data->date[dataIndex][0];
    counter = startDay-1;

    //Fill out each of the arrays
    for (int i=0; i<nDays; i++)
    {
        
        //If the years hasn't changed, increment counter
        if (data->date[dataIndex+i][0] == oldYear) counter++;
        //If it has changed, reset counter - this handles leap years
        else counter = 1;

        evap->day = counter;

        evap->P = asin(0.39795*cos(0.2163108 + 2.0 * atan(0.9671396*tan(0.00860*double(evap->day-186)))));
        evap->dayLength = 24.0 - (24.0/PI)*(acos((sin(0.8333*PI/180.0)+sin(data->gageLat*PI/180.0)*sin(evap->P))/(cos(data->gageLat*PI/180.0)*cos(evap->P))));
        evap->eStar = 0.6108*exp((17.27*data->avgTemp[dataIndex+i])/(237.3+data->avgTemp[dataIndex+i]));
        evap->PE[i] = (715.5*evap->dayLength*evap->eStar/24.0)/(data->avgTemp[dataIndex+i] + 273.2);

        oldYear = data->date[dataIndex+i][0];
    }

    return;
}