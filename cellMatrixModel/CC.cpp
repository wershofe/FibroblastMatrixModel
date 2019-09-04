// 19/10/16

#include <iostream>
#include "CC.h"
#include<ctime>
#include<cmath>
#include<cstdlib>
#include <stdio.h>
#include "Contact.h"
#include <vector>
#include<string>

using namespace std;

CC::CC(double radiusCD)
{
	setRadiusC(radiusCD);
	setAngleCI(0); 
	setNumberCI(-1); // = -1 if cc hasn't yet had a CI
}

CC::~CC()
{
}

void CC::setRadiusC(double radiusCD)
{
	radiusC = radiusCD;
}

void CC::setFirstPosition(int lengthD, double radiusCD, double thetaD) // random first position and orientation for each cc
{
	double randx = ((double) rand() / (RAND_MAX));
	double randy = ((double) rand() / (RAND_MAX));
	double rand01 = ((double) rand() / (RAND_MAX));
		
	Node n0;
	n0.setx(randx * lengthD);
	n0.sety(randy * lengthD);
	n0.setNodeType("H");
	nodes.push_back(n0);
	
	originalPosition.push_back(randx*lengthD);
	originalPosition.push_back(randy*lengthD);
	
	theta = thetaD;
}

void CC::setFirstPosition(double xD, double yD, double radiusCD, double thetaD) // dictated first position for a cc 
{
	double rand01 = ((double) rand() / (RAND_MAX));
	double randAngle = (rand01 * 2 * M_PI) ;
		
	Node n0;
	n0.setx(xD);
	n0.sety(yD);
	n0.setNodeType("H");
	nodes.push_back(n0);
	
	originalPosition.push_back(xD);
	originalPosition.push_back(yD);
	
	theta = thetaD;
}

double CC::getRadiusC()
{
	return radiusC;
}

void CC::setAngleCI(double angleCID)
{
	angleCI = angleCID;
}
double CC::getAngleCI()
{
	return angleCI;
}
void CC::setNumber(int numberD)
{
	number = numberD;
}
int CC::getNumber()
{
	return number;
}
void CC::setNumberCI(int numberCID)
{
	numberCI = numberCID;
}
int CC::getNumberCI()
{
	return numberCI;
}
std::vector<Node> CC::getNodes()
{
	return nodes;
}

void CC::setBinaryContact(int contactD) // Can be useful as a test to check the ccs were recognising contacts with each other
{
	binaryContact = contactD;
}
int CC::getBinaryContact()
{
	return binaryContact;
}
double CC::getTheta()
{
	return theta;
}
