// 28/10/16
#include <iostream>
#include "Parameter.h"
#include<cmath>
#include<cstdlib>
#include <stdio.h>
#include <vector>
#include <map>
#include<string>
#include<algorithm>
#include <fstream>

using namespace std;

Parameter::Parameter()
{

}
Parameter::~Parameter()
{

}

///////////////////////////////////////////////////////////////////
// sets
//////////////////////////////////////////////////////////////////
void Parameter::setAllParameters(std::string fileNameD)
{
	ifstream parameterFile(fileNameD.c_str());
	string dummyName;
	string variable;
	while(parameterFile >> dummyName >> variable)
	{
		if(dummyName == "numberVoxBox")
		{
			numberVoxBox = std::atoi(variable.c_str());
		}
		else if(dummyName == "rho")
		{
			rho = std::atoi(variable.c_str());	
		}
		else if(dummyName == "rhoCC")
		{
			rhoCC = std::atoi(variable.c_str());	
		}
		else if(dummyName == "speedMean")
		{
			speedMean = std::atof(variable.c_str());	
		}
		else if(dummyName == "speedSD")
		{
			speedSD = std::atof(variable.c_str());	
		}
		else if(dummyName == "polarity")
		{
			polarity = std::atof(variable.c_str());
		}
		else if(dummyName == "proliferation")
		{
			proliferation = std::atof(variable.c_str());
		}
		else if(dummyName == "radiusC")
		{
			radiusC = std::atof(variable.c_str());
		}
		else if(dummyName == "radiusCC")
		{
			radiusCC = std::atof(variable.c_str());
		}
		else if(dummyName == "CILa")
		{
			CILa = std::atof(variable.c_str());
		}
		else if(dummyName == "CILb")
		{
			CILb = std::atof(variable.c_str());
		}
		else if(dummyName == "hetCILa")
		{
			hetCILa = std::atof(variable.c_str());
		}
		else if(dummyName == "hetCILb")
		{
			hetCILb = std::atof(variable.c_str());
		}
		else if(dummyName == "CILtype")
		{
			CILtype = variable;
		}
		else if(dummyName == "w2")
		{
			w2 = std::atof(variable.c_str());
		}
		else if(dummyName == "w3")
		{
			w3 = std::atof(variable.c_str());
		}
		else if(dummyName == "w4")
		{
			w4 = std::atof(variable.c_str());
		}
		else if(dummyName == "w5")
		{
			w5 = std::atof(variable.c_str());
		}
		else if(dummyName == "length")
		{
			length = std::atoi(variable.c_str());
		}
		else if(dummyName == "trajectoryFileName")
		{
			trajectoryFileName = variable;
		}
		else if(dummyName == "fibreGridLength")
		{
			fibreGridLength = std::atoi(variable.c_str());
		}
		else if(dummyName == "boxesToRecord")
		{
			boxesToRecord = std::atoi(variable.c_str());
		}
		else if(dummyName == "sourceX")
		{
			sourceX = std::atoi(variable.c_str());
		}
		else if(dummyName == "sourceY")
		{
			sourceY = std::atoi(variable.c_str());
		}
		else if(dummyName == "depRate")
		{
			depRate = std::atoi(variable.c_str());
		}
		else if(dummyName == "reRate")
		{
			reRate = std::atoi(variable.c_str());
		}
		else if(dummyName == "degRate")
		{
			
			degRate = std::atoi(variable.c_str());
		}
		else if(dummyName == "IMT")
		{
			IMT = variable;
		}
		else if(dummyName == "IMI")
		{
			
			IMI = std::atoi(variable.c_str());
		}
		else if(dummyName == "CILdataHB")
		{
			readInCILdataHB(variable);
		}
		else if(dummyName == "maxFrame")
		{
			maxFrame = std::atoi(variable.c_str());
			
		}
		else if(dummyName == "overlapAllowance")
		{
			overlapAllowance = std::atof(variable.c_str());
		}
		else if(dummyName == "cis")
		{
			cis = std::atof(variable.c_str());
		}
		else if(dummyName == "fakeRoseSD")
		{
			fakeRoseSD = std::atof(variable.c_str());
		}
		else if(dummyName == "fakeRoseMean")
		{
			fakeRoseMean = std::atof(variable.c_str());
		}
		else if(dummyName == "aspectRatio")
		{
			aspectRatio = std::atoi(variable.c_str());	
		}
		else if(dummyName == "noBins")
		{
			noBins = std::atoi(variable.c_str());
		}
		else if(dummyName == "headSize")
		{
			headSize = std::atof(variable.c_str());
		}
		else if(dummyName == "probCIL")
		{
			probCIL = std::atof(variable.c_str());
		}
		else if(dummyName == "woundSize")
		{
			woundSize = std::atoi(variable.c_str());	
		}
	}
	w1 = 1-w2;
	if(w1 <0)
	{
		cout << "ERROR WAFFLE w1 < 0" << endl;
		exit(9);
	}
	
	setNumberVoxBox();
	setFibreGridLength();
	if(boxesToRecord==-1)
	{
		boxesToRecord = rho;
	}
	cout << "boxesToRecord = " << boxesToRecord << endl;
}

void Parameter::readInCILdataHB(std::string CILdataHBnameD)
{
	double angle;
	ifstream trajectoryFile(CILdataHBnameD);
	string line;
	while(getline(trajectoryFile,line))
	{
		CILdataHB.push_back(atof(line.c_str()));
	}
	if(CILdataHB.size() ==0)
	{
		cout << "Reading nothing!" << endl;
		exit(-1);
	}
}
void Parameter::setSourceX(int sourceXD)
{
	sourceX = sourceXD;
}
void Parameter::setSourceY(int sourceYD)
{
	sourceY = sourceYD;
}
void Parameter::setNumberVoxBox()
{
	numberVoxBox = floor(length/(2*aspectRatio*radiusC));;
}
void Parameter::setFibreGridLength()
{
	int temp=0;
	if(fibreGridLength == -1)
	{
		temp = ceil(length/radiusC);
		fibreGridLength = temp;
	}
	
	cout << "fibreGridLength = " << fibreGridLength << endl;
}

void Parameter::setRho(int rhoD)
{
	rho = rhoD;
}
void Parameter::setSpeedMean(double speedMeanD)
{
	speedMean = speedMeanD;
}
void Parameter::setSpeedVar(double speedSDD)
{
	speedSD = speedSDD;
}
void Parameter::setPolarity(double polarityD)
{
	polarity = polarityD;
}
void Parameter::setProliferation(double proliferationD)
{
	proliferation = proliferationD;
}
void Parameter::setRadiusC(double radiusCD)
{
	radiusC = radiusCD;
}
void Parameter::setRadiusCC(double radiusCCD)
{
	radiusCC = radiusCCD;
}
void Parameter::setCILa(double CILaD)
{
	CILa = CILaD;
}
void Parameter::setHetCILa(double hetCILaD)
{
	hetCILa = hetCILaD;
}
void Parameter::setHetCILb(double hetCILbD)
{
	hetCILb = hetCILbD;
}
void Parameter::setW1(double w1D)
{
	w1 = w1D;
}
void Parameter::setW2(double w2D)
{
	w2 = w2D;
}
void Parameter::setW3(double w3D)
{
	w3 = w3D;
}
void Parameter::setW4(double w4D)
{
	w4 = w4D;
}
void Parameter::setW5(double w5D)
{
	w5 = w5D;
}
void Parameter::setCILb(double CILbD)
{
	CILb = CILbD;
}
void Parameter::setCILtype(std::string CILtypeD)
{
	CILtype = CILtypeD;
}
void Parameter::setLength(int lengthD)
{
	length = lengthD;
}
void Parameter::setTrajectoryFileName(std::string trajectoryFileNameD)
{
	trajectoryFileName = trajectoryFileNameD;
}

///////////////////////////////////////////////////////////////////
// gets
//////////////////////////////////////////////////////////////////
int Parameter::getMaxFrame()
{
	return maxFrame;
}
int Parameter::getNumberVoxBox()
{
	return numberVoxBox;
}
int Parameter::getRho()
{
	return rho;
}
int Parameter::getRhoCC()
{
	return rhoCC;
}
double Parameter::getSpeedMean()
{
	return speedMean;
}
double Parameter::getSpeedSD()
{
	return speedSD;
}
double Parameter::getPolarity()
{
	return polarity;
}
double Parameter::getProliferation()
{
	return proliferation;
}
double Parameter::getRadiusC()
{
	return radiusC;
}
double Parameter::getRadiusCC()
{
	return radiusCC;
}
double Parameter::getCILa()
{
	return CILa;
}
double Parameter::getHetCILa()
{
	return hetCILa;
}
double Parameter::getHetCILb()
{
	return hetCILb;
}
double Parameter::getW1()
{
	return w1;
}
double Parameter::getW2()
{
	return w2;
}
double Parameter::getW3()
{
	return w3;
}
double Parameter::getW4()
{
	return w4;
}
double Parameter::getW5()
{
	return w5;
}
double Parameter::getCILb()
{
	return CILb;
}
std::string Parameter::getCILtype()
{
	return CILtype;
}
int Parameter::getLength()
{
	return length;
}
std::string Parameter::getTrajectoryFileName()
{
	return trajectoryFileName;
}
int Parameter::getFibreGridLength()
{
	return fibreGridLength;
}
int Parameter::getBoxesToRecord()
{
	return boxesToRecord;
}
int Parameter::getSourceX()
{
	return sourceX;
}
int Parameter::getSourceY()
{
	return sourceY;
}
int Parameter::getDepRate()
{
	return depRate;
}
int Parameter::getReRate()
{
	return reRate;
}
int Parameter::getDegRate()
{
	return degRate;
}
string Parameter::getIMT()
{
	return IMT;
}
int Parameter::getIMI()
{
	return IMI;
}
std::vector<double> Parameter::getCILdataHB()
{
	return CILdataHB;
}
double Parameter::getOverlapAllowance()
{
	return overlapAllowance;
}
double Parameter::getCis()
{
	return cis;
}
double Parameter::getFakeRoseSD()
{
	return fakeRoseSD;
}
double Parameter::getFakeRoseMean()
{
	return fakeRoseMean;
}
int Parameter::getAspectRatio()
{
	return aspectRatio;
}
int Parameter::getNoBins()
{
	return noBins;
}
double Parameter::getHeadSize()
{
	return headSize;
}
double Parameter::getProbCIL()
{
	return probCIL;
}
int Parameter::getWoundSize()
{
	return woundSize;
}
