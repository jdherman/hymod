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

#if defined(HYMOD) || defined(HBVMOD)

#ifndef __mopexdata_h
#define __mopexdata_h

#include <string>

using namespace std;

//MOPEX data structure
struct MOPEXData
{
    string ID;
    double gageLat;
    double gageLong;
    double DA;

    //Starting and ending dates of data to read
    int *dateStart;
    int *dateEnd;

    int nDays;          //Number of days of data
    int **date;         //Date of data [year, month, day]
    double *precip;     //Mean areal precipitation (mm)
    double *evap;       //Climatic potential evaporation (mm)
    double *flow;       //Streamflow discharge (mm)
    double *maxTemp;    //Maximum air temperature (Celsius) (should be daily)
    double *minTemp;    //Minimum air temperature (Celsius) (should be daily)
    double *avgTemp;    //Average air temperature (Celsius) (should be daily)
    //double *peAdjust;  //PE adjustment factors for each month
};

void readMOPEXData(MOPEXData *data, string filename);

#endif

#endif //HYMOD


