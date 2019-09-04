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

ofstream FileCompTrajtxt; 

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
	
	string name = "File_" + paramList.getTrajectoryFileName() + "_" + stringRand + ".txt";
	
	
	string nameNewFile = stringRand + "_parametersToRead.txt";
	ifstream inFile("parametersToRead.txt");
	ofstream outFile(nameNewFile);
	outFile << inFile.rdbuf();
	
	FileCompTrajtxt.open(name.c_str()); 

//-------------------------------------------------------------------------------------------------------------

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
			
			
		a.updateTME(frame);	
	}	
		cout << stringRand << endl;
	 	FileCompTrajtxt.close();
	  
	  t2=clock();
	  double diff ((double)t2-(double)t1);
	  return 0;		
}


//-------------------------------------------------------------------------------------------------------------------------


