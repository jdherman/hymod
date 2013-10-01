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

HyMod hymod;
HyMod *hymodPtr = &hymod;

//These are settings for a 10 year calibration period
int dateStart[3] = {1961, 10, 1};
int dateEnd[3]   = {1972, 9, 30};
int dayStartIndex = 274; //10-1 is day 274 of the year
//Note: For the calibration period specification, you should include the warmup period here
//The objective values get calculated on the time series values following the first year of simulation
int nDays = 4017; //11 years total - 1 year of warmup, 10 years of calibration, and 2 leap years.
int startingIndex = 5023-1; //The starting index of the data array corresponding with the start date

int PeriodLength;
int *Period;

//This function gets everything set up for running HyMod within e-NSGAII
void init_HyMod(int argc, char** argv)
{
    //Create pointers to difference parts of the HyMod structure
    MOPEXData *data        = &hymodPtr->data;
    HamonEvap *evap        = &hymodPtr->evap;
    HyModPars *pars        = &hymodPtr->pars;
    SnowDDPars *snowPars   = &hymodPtr->snowPars;
    HyModInStates *inState = &hymodPtr->inState;
    HyModModel *model      = &hymodPtr->model;
    SnowDDModel *snowModel = &hymodPtr->snowModel;

    pars->Nq = 1; // number of quickflow reservoirs
    pars->Kv = 1.0; // vegetation parameter
    snowPars->useSnowDD = true;

    //Allocate the data date arrays
    data->dateStart    = new int[3];
    data->dateEnd      = new int[3];
    data->dateStart[0] = dateStart[0];
    data->dateStart[1] = dateStart[1];
    data->dateStart[2] = dateStart[2];
    data->dateEnd[0]   = dateEnd[0];
    data->dateEnd[1]   = dateEnd[1];
    data->dateEnd[2]   = dateEnd[2];
    data->nDays = nDays;

	PeriodLength = nDays;
    Period = new int[PeriodLength];

    for (int d=0; d<PeriodLength; d++)
    {
        Period[d] = startingIndex+d;
    }

    //Read in input data
    string inputFile = argv[2];
    readMOPEXData(data, inputFile);

    //Calculate the Hamon Potential Evaporation for the time series
    calculateHamonPE(data, startingIndex, nDays, evap, dayStartIndex);

    //Allocate the model state arrays used by the snow model
    snowModel->snow      = new double [nDays];
    snowModel->melt      = new double [nDays];
    snowModel->store     = new double [nDays];
    
    //Allocate the model output arrays
    model->effPrecip     = new double [nDays];
    model->XHuz = new double [nDays];
    model->XCuz = new double [nDays];
    model->Xs   = new double [nDays];
    model->AE   = new double [nDays];
    model->OV   = new double [nDays];
    model->Qq   = new double [nDays];
    model->Qs   = new double [nDays];
    model->Q    = new double [nDays];
    model->Xq   = new double* [nDays];
    for (int i=0; i<nDays; i++) model->Xq[i] = new double [pars->Nq];

    return;
}

//This is the function that gets called whenever an individual gets evaluated
void calc_HyMod(double* parameters)
{
    //Create pointers to difference parts of the HyMod structure
    MOPEXData *data        = &hymodPtr->data;
    HamonEvap *evap        = &hymodPtr->evap;
    HyModPars *pars        = &hymodPtr->pars;
    SnowDDPars *snowPars   = &hymodPtr->snowPars;
    HyModInStates *inState = &hymodPtr->inState;
    HyModModel *model      = &hymodPtr->model;
    SnowDDModel *snowModel = &hymodPtr->snowModel;

    //Reinitialize everything to zero
    for (int k=0; k<nDays; k++)
    {
        snowModel->snow[k]  = 0.0;
        snowModel->melt[k]  = 0.0;
        snowModel->store[k] = 0.0;
        model->effPrecip[k] = 0.0;
        model->XHuz[k]      = 0.0;
        model->AE[k] = 0.0;
        model->OV[k] = 0.0;
        model->Qq[k] = 0.0;
        model->Qs[k] = 0.0;
        model->Q[k]  = 0.0;

        //These ones can optionally be initialized
        model->XCuz[k] = 0.0;
        model->Xs[k] = 0.0;
        for (int m=0; m<pars->Nq; m++) model->Xq[k][m] = 0.0;
    }

    //Put the variables associated with the individual into the HyMod parameter structure
    pars->Huz   = parameters[7];
    pars->B     = parameters[6];
    pars->alpha = parameters[5];
    pars->Kq    = parameters[1]; //1.0 parameters[3]; //Kq is now specified in days
    pars->Ks    = parameters[0]; //1.0 parameters[4]; //Ks is now specified in days

    snowPars->DDF = parameters[2];
    snowPars->Tth = parameters[4];
    snowPars->Tb  = parameters[3];

    //Run Model for Simulation Period
    HyModMo(*data, *evap, Period, PeriodLength, *pars, *snowPars, *inState, model, snowModel);

    //Calculate objective functions of calibration period (assumes first year is warmup period)
    // model->Q, data->flow, data->precip, evap->PE, startingIndex, PeriodLength
    return;
}

void deleteAllStuff(MOPEXData *data, HyModModel *model, SnowDDModel *snowModel)
{
	delete[] data->dateStart;
    delete[] data->dateEnd;

    delete[] snowModel->snow;
    delete[] snowModel->melt;
    delete[] snowModel->store;
    
    delete[] model->effPrecip;
    delete[] model->XHuz;
    delete[] model->XCuz;
    delete[] model->Xs;
    delete[] model->AE;
    delete[] model->OV;
    delete[] model->Qq;
    delete[] model->Qs;
    delete[] model->Q;
    
    //for (int i=0; i<nDays; i++) delete[] model->Xq[i];
	//delete[] model->Xq;
}

int main(int argc, char **argv)
{
	ifstream in; //Input file stream
    ofstream output;
    
	string sampleFile = argv[3];
	int numSamples = atoi(argv[4]);
	string outputFile = argv[5];

    int nParams = 8;
    double* parameters = new double [nParams]; 

    //Intialize
    init_HyMod(argc, argv);

	in.open(sampleFile.c_str(), ios_base::in);
	if (!in)
	{
		cout << "Error: sampleFile could not be opened!" << endl;
		exit(1);
	}

    //Loop through samples in Sobol file
    for (int s=0; s<numSamples; s++)
    {
		if (s%100 == 0) cout << s << endl;

        //Read in the parameter values from Sobol
        for (int i=0; i<nParams; i++)
        {
			in >> parameters[i];	
        }
		
		in.ignore(1000,'\n');

        //Run this sample
        calc_HyMod(parameters);
/*        of_decade << setw(15) << decadeObj->rmse[i];
		of_decade << endl;
*/    }

    //Close the output files
    //of_decade.close();

	deleteAllStuff(&hymodPtr->data, &hymodPtr->model, &hymodPtr->snowModel);

    return 0;
}