//13/6/16
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <stdio.h>
#include "FibreBox.h"
#include <vector>
#include <algorithm>

using namespace std;

FibreBox::FibreBox(int xD, int yD, int boxD, int dominantIndexD, int dominantDensityD, int depRateD, int reRateD, int degRateD, int noBinsD, bool woundedD)
{
	depositionRate = depRateD;
	rearrangementRate = reRateD;
	degradationRate = degRateD;
	noBins = noBinsD;
	wounded = woundedD;
	
	x = xD;
	y = yD;
	dominantIndex = dominantIndexD;
	dominantDensity = dominantDensityD;
	box = boxD;
	
	for(int i=0; i!=noBins; i++)
	{
		bins.push_back( i*(M_PI/noBins) );
		densityBins.push_back(0);
	}	
	densityBins[dominantIndex] = dominantDensity;
	dominantOrientation = bins[dominantIndexD];
}

FibreBox::FibreBox(int xD, int yD, int boxD, vector<int> densities, int depRateD, int reRateD, int degRateD, int noBinsD, bool woundedD)
{
	depositionRate = depRateD;
	rearrangementRate = reRateD;
	degradationRate = degRateD;
	noBins = noBinsD;
	wounded = woundedD;
	
	x = xD;
	y = yD;
	box = boxD;
	bins.reserve(densities.size());
	for(int i=0; i!=densities.size(); i++)
	{
		bins.push_back( i*(M_PI/noBins) );
		densityBins.push_back(densities[i]);
	}
	
	setDominantIndex();
	/*for(int i=0; i!=densities.size(); i++)
	{
		cout << bins[i] << "\t" <<densities[i] << endl;
	}
	cout << endl;*/	
}
FibreBox::~FibreBox()
{
}
void FibreBox::setX(int xD)
{
	x = xD;
}
void FibreBox::setY(int yD)
{
	y = yD;
}
void FibreBox::setBox(int boxD)
{
	box = boxD;
}

int FibreBox::getX()
{
	return x;
}
int FibreBox::getY()
{
	return y;
}
int FibreBox::getBox()
{
	return box;
}

//--------------------------------------------------

void FibreBox::setBins(std::vector<double> binsD)
{

}
void FibreBox::setDensityBins(std::vector<int> densityBinsD)
{

}
std::vector<double> FibreBox::getBins()
{
	return bins;
}
std::vector<int> FibreBox::getDensityBins()
{
	return densityBins;
}			
//---------------------------------------------------			
void FibreBox::processNewFibres(double thetaD)
{
	int index = computeBin(thetaD);
	degradeAll();
	rearrangeNeighbours(index);
	depositFibre(index);
	setDominantIndex();
}
int FibreBox::computeBin(double thetaD)
{
	int index;
	while(thetaD >M_PI)
	{
		thetaD -= M_PI;
	}
	vector<double> tempVec = bins;
	tempVec.push_back(M_PI);
	int tempIndex=0;
	double tempDistance=10000;
	for(int i=0; i!=tempVec.size(); i++)
	{	
		double newDist = (thetaD-tempVec[i])*(thetaD-tempVec[i]);
		if(newDist<tempDistance)
		{
			tempDistance=newDist;
			tempIndex=i;
		}
		else
		{
			break;
		}
	}
	if(tempIndex==(tempVec.size()-1))
	{
		tempIndex=0;
	}
	index=tempIndex;
	return index;
}
void FibreBox::degradeAll()
{
	for(int i=0; i!=densityBins.size(); i++)
	{
		degradeBin(i);
	}
}
void FibreBox::degradeBin(int i)
{
	if(densityBins[i]>0)
	{
		densityBins[i]-=degradationRate;
	}
	if(densityBins[i]<0)
	{
		densityBins[i]=0;
	}
}
void FibreBox::rearrangeBin(int i)
{
	if(densityBins[i]>0)
	{
		densityBins[i]-=rearrangementRate;
	}
	if(densityBins[i]<0)
	{
		densityBins[i]=0;
	}
}
void FibreBox::rearrangeNeighbours(int indexD)
{
	int neighbourLeft = indexD-1;
	int neighbourRight = indexD+1;
	
	
	if(neighbourLeft<0)
	{
		neighbourLeft+=bins.size();
	}
	if(neighbourRight>=bins.size())
	{
		neighbourRight-=bins.size();
	}
	
	int totalRearranged=0;
	
	if(densityBins[neighbourLeft]>= rearrangementRate)
	{
		totalRearranged+=rearrangementRate;
	}
	else
	{
		totalRearranged+=densityBins[neighbourLeft];
	}
	
	if(densityBins[neighbourRight]>= rearrangementRate)
	{
		totalRearranged+=rearrangementRate;
	}
	else
	{
		totalRearranged+=densityBins[neighbourRight];
	}
	densityBins[indexD]+=totalRearranged;
	
	rearrangeBin(neighbourLeft);
	rearrangeBin(neighbourRight);
	
}
void FibreBox::depositFibre(int indexD)
{
	densityBins[indexD] += depositionRate; // for now depositing 10 fibres in direction of CAF
}
void FibreBox::setDominantIndex()
{	
	int totalDensity=0;
	
	for(int i=0; i!=densityBins.size(); i++)
	{
		totalDensity +=densityBins[i];
	}
	double randNo = ((double) rand() / (RAND_MAX))*totalDensity;
	
	int tempIndex=0;
	int scroll=densityBins[tempIndex];
	
	while(scroll<randNo)
	{
		tempIndex++;
		scroll += densityBins[tempIndex];
	}
	if(totalDensity==0)
	{
		int randNo = ((double) rand() / (RAND_MAX))*densityBins.size();
		tempIndex = randNo;
	}
	
	dominantIndex=tempIndex;
	dominantOrientation=bins[dominantIndex];
	dominantDensity=densityBins[dominantIndex];
}
int FibreBox::getDominantIndex()
{
	return dominantIndex;
}
void FibreBox::setDominantDensity(int dominantDensityD)
{
	dominantDensity = dominantDensityD;
}
void FibreBox::setDominantOrientation(double dominantOrientationD)
{
	dominantOrientation = dominantOrientationD;
}
int FibreBox::getDominantDensity()
{
	return dominantDensity;
}
double FibreBox::getDominantOrientation()
{
	return dominantOrientation;
}
bool FibreBox::getWounded()
{
	return wounded;
}
