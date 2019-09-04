// 4/11/16
#include "Fibroblast.h"
#include <vector>
#include <cmath>
#include "TME.h"
#include "Contact.h"
#include <stdio.h>
#include <stdarg.h>
#include <vector>
#include <math.h>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string.h>
#include <string>
#include <ctime>
#include <map>
#include "Node.h"
#include "Parameter.h"
#include "FibreBox.h"
#include "CC.h" 
#include "CAT.hpp"
#include <omp.h>

using namespace std;

// --------------------------------------------------------------------------------------------------------------
ofstream FileCompFibreBoxtxt;
ofstream FileCompCC;
ofstream FileCAT;
ofstream FileCompTrajtxt; 
ofstream FileGreyScaletxt;
ofstream FileQuantFibrestxt;
ofstream FileNewMatrix;
ofstream FileOldMatrix;

Parameter paramList;
int frame;
int maxFrame;
TME a;  // a is the key to access all thgridBoxLength 150
void makeSourceSink(vector<iPair> &source, vector<iPair> &sink, int xc, int yc, int radius);
void plot(double chemoattr[][nblock], double t_now);
// --------------------------------------------------------------------------------------------------------------
 
int main(int argc, char* argv[])
{
	clock_t t1,t2;
	t1=clock();

	int i;
	ifstream random("/dev/urandom",ios::binary); // doing a better randomness - Raph's code
	unsigned int seed;
	seed=(unsigned int)random.get();
	for(i=0; i<sizeof(unsigned int)-1;i++)
	{
		seed=seed<<8;
		seed|=(unsigned int)random.get();
	}
	random.close();
	srand(seed);

//------------------------------------------------------------------------------------------------------------------	
	paramList.setAllParameters("parametersToRead.txt");
	maxFrame=paramList.getMaxFrame();
	/*if(paramList.getProliferation()!=0)
	{
		paramList.setRho(50);
	}
	else
	{
		//paramList.setRho(500); removed for now - put back in soon
	}*/
	a.readParameters(paramList);
	string stringRand = static_cast<ostringstream*>( &(ostringstream() << seed) )->str();
	cout << stringRand << endl;
	string stringPolarity = static_cast<ostringstream*>( &(ostringstream() << paramList.getPolarity()) )->str();
	string stringW2 = static_cast<ostringstream*>( &(ostringstream() << paramList.getW2()) )->str();
	string stringW5 = static_cast<ostringstream*>( &(ostringstream() << paramList.getW5()) )->str();
	string stringDepRate = static_cast<ostringstream*>( &(ostringstream() << paramList.getDepRate()) )->str();
	string stringReRate = static_cast<ostringstream*>( &(ostringstream() << paramList.getReRate()) )->str();
	string stringDegRate = static_cast<ostringstream*>( &(ostringstream() << paramList.getDegRate()) )->str();
	string name = "File_" + paramList.getTrajectoryFileName() + "_" + stringRand + ".txt";
	string nameGreyScale = name + "GreyScale.txt";
	string nameOldMatrix = stringRand + "_OldMatrix.txt";
	
	string nameFibreBoxtxt = name + "_FibreBox.txt";
	string nameQuantFibrestxt = name + "QuantFibres.txt";
	string nameNewFile = stringRand + "_parametersToRead.txt";
	string nameNewMatrix = stringRand + "_NewMatrix.txt";
	ifstream inFile("parametersToRead.txt");
	ofstream outFile(nameNewFile);
	outFile << inFile.rdbuf();
	
	FileCompFibreBoxtxt.open(nameFibreBoxtxt.c_str());
	FileQuantFibrestxt.open(nameQuantFibrestxt.c_str()); 
	FileNewMatrix.open(nameNewMatrix.c_str());
	FileOldMatrix.open(nameOldMatrix.c_str());
	FileCompTrajtxt.open(name.c_str()); 
	FileGreyScaletxt.open(nameGreyScale.c_str());
//-------------------------------------------------------------------------------------------------------------
	 double xdim = 1., ydim = 1. ;
    
	    vector<iPair> source, sink;
	    int xc = nblock/2, yc = nblock/2, radius = nblock/10;
	    //source.push_back({2*nblock/5, nblock/2});
	    //source.push_back({3*nblock/5, nblock/2});
	    makeSourceSink(source, sink, xc, yc, radius);
	  //  xc = 3*nblock/5;
	   // makeSourceSink(source, sink, xc, yc, radius);
	    
	    double chemoattr[nblock][nblock] = {};
	    
	    double t = 0.;
	    
	    int nSaveFig = N_DT_SAVE_FIG;
	    int cnt = 0;
	    while (t<=tt)
	    {
		singleStep(chemoattr, xdim, ydim, source, sink);
		
		if (cnt % nSaveFig == 0)
		{
		    cout << "Making CAT. time now: t = " << t << endl;
		    plot(chemoattr, t);
		}
		
		t += dt;
		cnt++;
	    }
//-------------------------------------------------------------------------------------------------------------
	frame = 0;
	a.initialiseTME(); 
	while(frame !=maxFrame)
	{	
		cout << "\t Frame " << frame << endl;	
		//if(frame>0)
		//{			 
			vector<Fibroblast> p = a.getFibroblasts();
			int noCells = p.size();
			
			//if(frame%700==0)
			//if(frame==maxFrame-1)
			{
				for(int i=0; i!= p.size(); i++)
				{
					int number = p[i].getNumber();
					double x = p[i].getNodes()[0].getx();
					double y = p[i].getNodes()[0].gety();
					double theta = p[i].getTheta();
					int mobility = p[i].getMobility();
					double CI = p[i].getAngleCI();
				
					FileCompTrajtxt << number << "\t" << frame << "\t" << x << "\t"
					 << y << "\t" << theta << endl;	
				}
			}
			if(frame==maxFrame-1)
			{
				for(int i=0; i!= p.size(); i++)
				{
					int number = p[i].getNumber();
					double x = p[i].getNodes()[0].getx();
					double y = p[i].getNodes()[0].gety();
					double theta = p[i].getTheta();
					int contacts = p[i].getCurrentContacts().size();
				
					FileGreyScaletxt << number << "\t" << frame << "\t" << x << "\t"
					 << y << "\t" << theta << "\t" << contacts << endl;
				 }
			}
			
			//if(frame%700==0)
			if(frame==maxFrame-1)
			{
				int boxesToRecord = paramList.getBoxesToRecord();

				vector<FibreBox> quantFib = a.getQuantBoxes();
				for(int i=0; i!=quantFib.size(); i++)
				{
					int box = quantFib[i].getBox();
					int x = quantFib[i].getX(); 
					int y = quantFib[i].getY(); 
					double orientation = quantFib[i].getDominantOrientation();
					int density = quantFib[i].getDominantDensity();
					// These two are just so that can process ECM same as Cells - this is rubbish code, putting in dummy values
					int mobility = 0;
					double CI = 0;
				
					double xCoord = x*((double)paramList.getLength()/paramList.getFibreGridLength());
					double yCoord = y*((double)paramList.getLength()/paramList.getFibreGridLength());
				
					FileQuantFibrestxt << box << "\t" << frame << "\t" << xCoord << "\t"
					 << yCoord << "\t" << orientation << endl;
				
				}
			}			
			
			//if(frame%700==0)
			if(frame==maxFrame-1)
			{
				vector<FibreBox> fib = a.getFibreBoxes();
				
				int noBoxes = fib.size();

				for(int i=0; i!= fib.size(); i++)
				{	
					int x = fib[i].getX(); 
					int y = fib[i].getY(); 
					int box = fib[i].getBox();
					int density = fib[i].getDominantDensity();
					double orientation = fib[i].getDominantOrientation(); 
					
					FileCompFibreBoxtxt << box << "\t" << frame << "\t" << x << "\t" << y << "\t" << orientation << "\t" << density << endl;
				} 
				FileCompFibreBoxtxt << endl;
			}
			
			if(frame%700==0)
			//if(frame==maxFrame-1)
			{
				vector<FibreBox> fib = a.getFibreBoxes();
				int noBoxes = fib.size();
			
				for(int i=0; i!= fib.size(); i++)
				{	
					
						int x = fib[i].getX(); 
						int y = fib[i].getY(); 
						int box = fib[i].getBox();
						int density =0;
						if(fib[i].getWounded()==true)
						{
							density = fib[i].getDominantDensity();
						}
					
						double orientation = fib[i].getDominantOrientation(); 
						
						FileNewMatrix << box << "\t" << frame << "\t" << x << "\t" << y << "\t" << orientation << "\t" << density << endl;
					
				} 
			}
			
			if(frame==maxFrame-1)
			{
				vector<FibreBox> fib = a.getFibreBoxes();
				int noBoxes = fib.size();
			
				for(int i=0; i!= fib.size(); i++)
				{	
					int x = fib[i].getX(); 
					int y = fib[i].getY(); 
					int box = fib[i].getBox();
					std::vector<int> densityBins = fib[i].getDensityBins();
					FileOldMatrix << box << "\t" << x << "\t" << y << "\t" ;
					for(int j=0; j!=densityBins.size(); j++)
					{
						FileOldMatrix << densityBins[j] << "\t";
					}
					
					FileOldMatrix << endl;
				} 
				FileOldMatrix << endl;
			}
			
		
		/*
		vector<CC> CC = a.getCCs();
		int noCCs = CC.size();
		FileCompCC.write(reinterpret_cast<char*>(&noCCs),sizeof(int));
		for(int i=0; i!= CC.size(); i++)
		{	
			int number = CC[i].getNumber();
			double x = CC[i].getNodes()[0].getx();
			double y = CC[i].getNodes()[0].gety(); 	
			
			FileCompCC.write(reinterpret_cast<char*>(&number),sizeof(int));
			FileCompCC.write(reinterpret_cast<char*>(&frame),sizeof(int));
			FileCompCC.write(reinterpret_cast<char*>(&x),sizeof(double));
			FileCompCC.write(reinterpret_cast<char*>(&y),sizeof(double));
		} */
		/*
		vector<double> Mij = a.getMij();
		for(int i=0; i!=Mij.size(); i++)
		{
			double cat = Mij[i];
			
			FileCAT.write(reinterpret_cast<char*>(&i),sizeof(int)); 
			FileCAT.write(reinterpret_cast<char*>(&frame),sizeof(int));
			FileCAT.write(reinterpret_cast<char*>(&cat),sizeof(double));
		}
		*/
		//}
		a.updateTME(frame);	
	}	
	cout << stringRand << endl;
	  FileCompFibreBoxtxt.close();
	  FileQuantFibrestxt.close();
	  FileCompCC.close();
	  FileCompTrajtxt.close();
	  FileGreyScaletxt.close();
	  FileNewMatrix.close();
	  FileOldMatrix.close();
	 // FileCAT.close();
	  t2=clock();
	  double diff ((double)t2-(double)t1);
	  return 0;		
}


//-------------------------------------------------------------------------------------------------------------------------

void makeSourceSink(vector<iPair> &source, vector<iPair> &sink, int xc, int yc, int radius)
{
    // (1) circular source pixels
    cout << "making source ... " << endl;
    //int radius = nblock / 20;
    //int xc = nblock/2, yc = nblock/2;
    int xx, xmax = xc + radius;
    int yy = yc - radius, ymax = yc + radius;
    while (yy++ < ymax)
    {
        xx = xc - radius;
        while (xx++ < xmax)
        {
            double d2c = sqrt( (xx-xc)*(xx-xc) + (yy-yc)*(yy-yc) );
            if (d2c < radius)
            {
                source.push_back({xx,yy});
                cout << xx << " " << yy << endl;
            }
            
        }
    }
    
    // (2) randomly located sink pixels
    cout << "making sink ... " << endl;
    yy = 0;
    xmax = nblock-1;
    ymax = nblock-1;
    while (yy++ < ymax)
    {
        xx = 0;
        while (xx++ < xmax)
        {
            double d2c = sqrt( (xx-xc)*(xx-xc) + (yy-yc)*(yy-yc) );
            if (d2c > radius*2)
            {
                double r = ((double) rand() / (RAND_MAX));
                if (r > 0.9999)
                {
                    sink.push_back({xx,yy});
                    cout << xx << " " << yy << endl;
                }
                
            }
        }
    }
}


void plot(double chemoattr[][nblock], double t_now)
{
    // -------------- write to a temp .dat file --------------
    ofstream ca_file;
    string nameTime = static_cast<ostringstream*>( &(ostringstream() << t_now) )->str();
    string nameCAT = "mCAT/mCAT_" + nameTime + ".txt"; 
    ca_file.open(nameCAT.c_str());
    for (int i = 0; i < nblock; i ++)
    {
        for (int j = 0; j < nblock; j ++)
        {
            ca_file << chemoattr[i][j] << "\t";
        }
        ca_file << "\n";
    }
    ca_file.close();    
}

