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

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <sstream>
#include <algorithm>
#include <time.h>
#include <iomanip>
#include <math.h>

#include "MOPEXData.h"

using namespace std;
const double PI = 3.141592653589793238462;

struct hymod_parameters
{
    //User specified parameters
    double Huz;      //Maximum height of soil moisture accounting tank - Range [0, Inf]
    double B;        //Scaled distribution function shape parameter    - Range [0, 2]
    double alpha;    //Quick/slow split parameter                      - Range [0, 1]
    int    Nq;       //Number of quickflow routing tanks               - Range [1, Inf] (typically set to <3)
    double Kq;       //Quickflow routing tanks' rate parameter         - Range [0, 1]
    double Ks;       //Slowflow routing tank's rate parameter          - Range [0, 1]

    // Snow parameters (degree-day model)
    double DDF;        //Degree day factor                        - Range [ 0, 2]
    double Tth;        //Temperature threshold                    - Range [-5, 5] 
    double Tb;         //Base temperature to calculate melt       - Range [-5, 5]

    // Given/calculated parameters
    double Kv;       //Vegetation adjustment to PE                     - Range [0, 2]
    double Cpar;     //Maximum combined contents of all stores (calculated from Huz and b)
};

struct hymod_states
{
    double *XHuz;        //Model computed upper zone soil moisture tank state height
    double *XCuz;        //Model computed upper zone soil moisture tank state contents
    double **Xq;         //Model computed quickflow tank states contents
    double *Xs;          //Model computed slowflow tank state contents
    double *snow_store;      //State of snow reservoir
};

struct hymod_fluxes
{
    double *effPrecip;   //Effective rain entering the SMA model (melt+precip if using snow model)
    double *AE;          //Model computed actual evapotranspiration flux
    double *OV;          //Model computed precipitation excess flux
    double *Qq;          //Model computed quickflow flux
    double *Qs;          //Model computed slowflow flux
    double *Q;           //Model computed total streamflow flux
    double *snow;        //Daily snow
    double *melt;        //Snow melt
    double *PE;          // Potential ET (Hamon)
};

struct HyMod
{
    MOPEXData data;
    hymod_parameters parameters;
    hymod_states states;   
    hymod_fluxes fluxes;
};

// Global hymod structure, accessible by all files
extern struct HyMod hymod;

//Function Prototypes
void zero_states_and_fluxes(int ndays);
void hymod_allocate(int ndays);
void hymod_delete(int ndays);
void PDM_soil_moisture(int modelDay, int dataDay);
double Nash(double K, int N, double Qin, double *X);
double snowDD(int modelDay, int dataDay);
void calculateHamonPE(int dataIndex, int nDays, int startDay);
