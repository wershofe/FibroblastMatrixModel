//9/5/16
#include <iostream>
#include<ctime>
#include<cmath>
#include<cstdlib>
#include <stdio.h>
#include "Contact.h"


using namespace std;

Contact::Contact()
{
}

Contact::~Contact()
{
}
////////////////////////////////////////////////////////////////////
// SETS //
///////////////////////////////////////////////////////////////////
void Contact::setNodeTypeI(std::string nodeTypeD)
{
	nodeTypeI = nodeTypeD;
}
void Contact::setNodeTypeJ(std::string nodeTypeD)
{
	nodeTypeJ = nodeTypeD;
}
void Contact::setAngleI(double angleID)
{
	angleI = angleID;
}
void Contact::setAngleJ(double angleJD)
{
	angleJ = angleJD;
}
void Contact::setAngleDifference()
{
	angleDifference = abs(angleI - angleJ);
	while(angleDifference > M_PI)
	{
		angleDifference = (2*M_PI) - angleDifference; // So that difference angle is never more than pi
	}
}

void Contact::setCellI(int CellD)
{
	cellI = CellD;
}
void Contact::setCellJ(int CellD)
{
	cellJ = CellD;
}

////////////////////////////////////////////////////////////////////
// GETS //
///////////////////////////////////////////////////////////////////
std::string Contact::getNodeTypeI()
{
	return nodeTypeI;
}
std::string Contact::getNodeTypeJ()
{
	return nodeTypeJ;
}
double Contact::getAngleI()
{
	return angleI;
}
double Contact::getAngleJ()
{
	return angleJ;
}
double Contact::getAngleDifference()
{
	return angleDifference;
}

int Contact::getCellI()
{
	return cellI;
}
int Contact::getCellJ()
{
	return cellJ;
}












