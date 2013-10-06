/*
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
*/

#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>

using namespace std;

struct MOPEXData
{
    string ID;
    double gageLat;
    double gageLong;
    double DA;

    int nDays;          //Number of days of data
    int **date;         //Date of data [year, month, day]
    double *precip;     //Mean areal precipitation (mm)
    double *evap;       //Climatic potential evaporation (mm)
    double *flow;       //Streamflow discharge (mm)
    double *maxTemp;    //Maximum air temperature (Celsius) (should be daily)
    double *minTemp;    //Minimum air temperature (Celsius) (should be daily)
    double *avgTemp;    //Average air temperature (Celsius) (should be daily)
};

//Function to read in the MOPEX data (precip, flow, temp, AE, etc.)
void readMOPEXData(MOPEXData *data, string filename);
