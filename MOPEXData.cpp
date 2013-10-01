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
#include <string>
#include <fstream>
#include <cstdlib>

#include "MOPEXData.h"

using namespace std;

//Function to read in the MOPEX data (precip, flow, temp, AE, etc.)
void readMOPEXData(MOPEXData *data, string filename)
{

    ifstream in;
    string sJunk = "";
    int ijunk;
    double dTemp;

    in.open(filename.c_str(), ios_base::in);
    if(!in)
    {
        cout << "The input file specified: " << filename << " could not be found!" << endl;
        exit(1);
    }

    //Look for the <GAGE_ID> key
    while (sJunk != "<GAGE_ID>")
    {
        in >> sJunk;
    }
    in >> data->ID;
    //Return to the beginning of the file
    in.seekg(0, ios::beg);

    //Look for the <GAGE_LATITUDE> key
    while (sJunk != "<GAGE_LATITUDE>")
    {
        in >> sJunk;
    }
    in >> data->gageLat;
    //Return to the beginning of the file
    in.seekg(0, ios::beg);

    //Look for the <GAGE_LONGITUDE> key
    while (sJunk != "<GAGE_LONGITUDE>")
    {
        in >> sJunk;
    }
    in >> data->gageLong;
    //Return to the beginning of the file
    in.seekg(0, ios::beg);

    //Look for the <DRAINAGE_AREA> key
    while (sJunk != "<DRAINAGE_AREA>")
    {
        in >> sJunk;
    }
    in >> data->DA;
    //Return to the beginning of the file
    in.seekg(0, ios::beg);

    //Look for the <TIME_STEPS> key
    while (sJunk != "<TIME_STEPS>")
    {
        in >> sJunk;
    }
    in >> data->nDays;
    //Return to the beginning of the file
    in.seekg(0, ios::beg);

    //Allocate the arrays
    data->date = new int* [data->nDays];
    for (int i=0; i<data->nDays; i++) data->date[i] = new int[3];
    data->precip   = new double[data->nDays];
    data->evap     = new double[data->nDays];
    data->flow     = new double[data->nDays];
    data->maxTemp  = new double[data->nDays];
    data->minTemp  = new double[data->nDays];
    data->avgTemp  = new double[data->nDays];

    //Experimenting with reading in PE adjustment factors - Not tested
    //data->peAdjust = new double[12];

    ////Look for the <PE_ADJUST_START> key
    //while (sJunk != "<PE_ADJUST_START>")
    //{
    //    in >>sJunk;
    //}
    //in.ignore(1000,'\n');
    //for (int i=0; i<12; i++)
    //{
    //    in >> ijunk >> data->peAdjust[i];
    //    in.ignore(1000,'\n');
    //}
    ////Return to the beginning of the file
    //in.seekg(0, ios::beg);

    //Look for the <DATA_START> key
    while (sJunk != "<DATA_START>")
    {
        in >> sJunk;
    }
    
    //Once we found the key, ignore the rest of the line and move to the data
    in.ignore(1000,'\n');

    //Loop through all of the input data and read in this order:
    for (int i=0; i<data->nDays; i++)
    {
        in >> dTemp;
        data->date[i][0] = int(dTemp);
        in >> dTemp;
        data->date[i][1] = int(dTemp);
        in >> dTemp;
        data->date[i][2] = int(dTemp);
        in >> data->precip[i] >> data->evap[i] >> data->flow[i] >> data->maxTemp[i] >> data->minTemp[i];
        in.ignore(1000,'\n');
        //While we're at it, calculate average T
        data->avgTemp[i] = (data->maxTemp[i] + data->minTemp[i])/2.0;
    }

    //Close the input file
    in.close();

    return;
}
