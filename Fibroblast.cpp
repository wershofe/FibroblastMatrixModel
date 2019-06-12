  // 19/10/16

#include <iostream>
#include "Fibroblast.h"
#include<ctime>
#include<cmath>
#include<cstdlib>
#include <stdio.h>
#include "Contact.h"
#include <vector>
#include<string>

using namespace std;

Fibroblast::Fibroblast(double radiusCD, double speedMeanD, double speedSDD)
{
	setRadiusC(radiusCD);
	setSpeed(speedMeanD, speedSDD);
	setAngleCI(0); 
	setNumberCI(-1); // = -1 if fibroblast hasn't yet had a CI
	setAttractionAngle(0); // dummy value so don't work with uninitialised value in determineNewAngle()
	setAttractionAmplitude(0); // dummy value so don't work with uninitialised value in determineNewAngle()
	setTempNewAngle(0.0);
	setMobility(1);
}

Fibroblast::~Fibroblast()
{
}

void Fibroblast::setRadiusC(double radiusCD)
{
	radiusC = radiusCD;
}

void Fibroblast::setFirstPosition(int lengthD, double radiusCD, int aspectRatioD) // random first position and orientation for each fibroblast
{
	nodes.clear();
	overUnder = "nothing";
	nucBarrier = -1;
	double randx = ((double) rand() / (RAND_MAX));
	double randy = ((double) rand() / (RAND_MAX));
	double rand01 = ((double) rand() / (RAND_MAX));
	double randAngle = (rand01 * 2 * M_PI) ;
	theta = randAngle;
	double c = cos(theta);
	double s = sin(theta);
	
	Node n0;
	n0.setx(randx * lengthD);
	n0.sety(randy * lengthD);
	n0.setNodeType("H");
	nodes.push_back(n0);
	
	originalPosition.push_back(randx*lengthD);
	originalPosition.push_back(randy*lengthD);

	for(int i=0; i!=aspectRatioD; i++)
	{
		Node n1;
		n1.setx((randx * lengthD)- ((2*i)+1.5*radiusCD*c)); // 1.5 account for the head being half the size
		n1.sety((randy * lengthD)- ((2*i)+1.5*radiusCD*s));
		n1.setNodeType("B");
		nodes.push_back(n1);
	}
	double j=1.5 + 2*(aspectRatioD-1) + 1.5;
	Node n2;
	n2.setx((randx * lengthD)- j*radiusCD*c);
	n2.sety((randy * lengthD)- j*radiusCD*s);
	n2.setNodeType("T");
	nodes.push_back(n2);
	
	Node n4; // this is the "Centre" of the fibroblast - not graphically shown, about which the other nodes are positioned
	n4.setx((randx * lengthD)- (j/2)*radiusCD*c);
	n4.sety((randy * lengthD)- (j/2)*radiusCD*s);
	n4.setNodeType("C");
	nodes.push_back(n4);	
}

void Fibroblast::setFirstPosition(double xD, double yD, double radiusCD, int aspectRatioD) // dictated first position for a fibroblast 
{
	nodes.clear();
	overUnder = "nothing";
	nucBarrier = -1;
	double rand01 = ((double) rand() / (RAND_MAX));
	double randAngle = (rand01 * 2 * M_PI) ;
	theta = randAngle;
	double c = cos(theta);
	double s = sin(theta);
	
	Node n0;
	n0.setx(xD);
	n0.sety(yD);
	n0.setNodeType("H");
	nodes.push_back(n0);
	
	originalPosition.push_back(xD);
	originalPosition.push_back(yD);

	for(int i=0; i!=aspectRatioD; i++)
	{
		Node n1;
		n1.setx(xD- ((2*i)+1.5*radiusCD*c)); // 1.5 account for the head being half the size
		n1.sety(yD- ((2*i)+1.5*radiusCD*s));
		n1.setNodeType("B");
		nodes.push_back(n1);
	}
	
	double j=1.5 + 2*(aspectRatioD-1) + 1.5;
	Node n2;
	n2.setx(xD- j*radiusCD*c);
	n2.sety(yD- j*radiusCD*s);
	n2.setNodeType("T");
	nodes.push_back(n2);
	
	Node n4; // this is the "Centre" of the fibroblast - not graphically shown, about which the other nodes are positioned
	n4.setx(xD- (j/2)*radiusCD*c);
	n4.sety(yD- (j/2)*radiusCD*s);
	n4.setNodeType("C");
	nodes.push_back(n4);	
}

void Fibroblast::setFirstPosition(double xD, double yD, double thetaD, double radiusCD, int aspectRatioD)// dictated first position and orientation for a fibroblast 
{
	
	nodes.clear();
	overUnder = "nothing";
	nucBarrier = -1;
	
	theta = thetaD;
	double c = cos(thetaD);
	double s = sin(thetaD);

	Node n0;
	n0.setx(xD);
	n0.sety(yD);
	n0.setNodeType("H");
	nodes.push_back(n0);
	
	originalPosition.push_back(xD);
	originalPosition.push_back(yD);

	for(int i=0; i!=aspectRatioD; i++)
	{
		Node n1;
		n1.setx(xD- ((2*i)+1.5*radiusCD*c)); // 1.5 account for the head being half the size
		n1.sety(yD- ((2*i)+1.5*radiusCD*s));
		n1.setNodeType("B");
		nodes.push_back(n1);
	}
	
	double j=1.5 + 2*(aspectRatioD-1) + 1.5;
	Node n2;
	n2.setx(xD- j*radiusCD*c);
	n2.sety(yD- j*radiusCD*s);
	n2.setNodeType("T");
	nodes.push_back(n2);
	
	Node n4; // this is the "Centre" of the fibroblast - not graphically shown, about which the other nodes are positioned
	n4.setx(xD- (j/2)*radiusCD*c);
	n4.sety(yD- (j/2)*radiusCD*s);
	n4.setNodeType("C");
	nodes.push_back(n4);	
	
	
}

double Fibroblast::getTheta()
{
	if(theta>2*M_PI)
	{
		theta-=2*M_PI;
	}
	if(theta<0)
	{
		theta+=2*M_PI;
	}
	return theta;
}

void Fibroblast::setSpeed(double speedD)
{
	speed = speedD;
}

void Fibroblast::setSpeed(double speedMeanD, double speedSDD)
{
	double u1 = (double) rand() / (RAND_MAX);
	double u2 = (double) rand() / (RAND_MAX);
	double z0 = sqrt(-2*log(u1))*cos(2*M_PI*u2); // box-muller transform
	speed = abs(speedSDD*z0 + speedMeanD);
}

void Fibroblast::setAngle(double angleD) // factoring in the noise of the system
{
	theta = angleD;
}

void Fibroblast::setPosition(double radiusCD, int lengthD, int aspectRatioD)
{
	const double c = cos(theta);
	const double s = sin(theta);
	
	double x = speed * c * mobility; //14-2-17 mobility added in.
	double y = speed * s * mobility;
	int anchor; // establish which node is the centre
	
	
	// set centre
	if(nodes[nodes.size()-1].getNodeType() != "C" || nodes[nodes.size()-2].getNodeType() != "T" || nodes[0].getNodeType() != "H")
	{
		cout << "ERROR SPAGHETTI!" << endl;
		exit(-30);
	}
	anchor = nodes.size()-1;
	nodes[anchor].setx(nodes[anchor].getx() + x);
	nodes[anchor].sety(nodes[anchor].gety() + y);

	//set head
	double h = aspectRatioD+0.5;
	
	nodes[0].setx(nodes[anchor].getx() + h*radiusCD*c);
	nodes[0].sety(nodes[anchor].gety() + h*radiusCD*s);
	
	// set tail
	nodes[nodes.size()-2].setx(nodes[anchor].getx() - h*radiusCD*c);
	nodes[nodes.size()-2].sety(nodes[anchor].gety() - h*radiusCD*s);
	
	// need to work out if AR is odd or even, set position of nodes accordingly
	bool even = false;
	if(aspectRatioD %2 ==0)
	{
		even = true;
	}
	if(even == true)
	{
		double j= aspectRatioD-1; // position of first body node
		for(int i=1; i!=(nodes.size()-2); i++)
		{			
			nodes[i].setx(nodes[anchor].getx() +j*radiusCD*c);
			nodes[i].sety(nodes[anchor].gety() +j*radiusCD*s);
			j-=2;
		}
	}
	else // if aspect ratio (ie number of body nodes) is odd
	{	
		double j= aspectRatioD-1; // position of first body node
		for(int i=1; i!=(nodes.size()-2); i++)
		{	
			nodes[i].setx(nodes[anchor].getx() +j*radiusCD*c);
			nodes[i].sety(nodes[anchor].gety() +j*radiusCD*s);
			j-=2;
		}
	}
}

void Fibroblast::applyBoundaryConditions(int lengthD)
{
	for(int i=0; i!=nodes.size(); i++)
	{
		while(nodes[i].getx() < 0){nodes[i].setx(nodes[i].getx() + lengthD);} //Left
		while(nodes[i].getx() > lengthD){nodes[i].setx(nodes[i].getx() - lengthD);} //Right
		while(nodes[i].gety() < 0){nodes[i].sety(nodes[i].gety() + lengthD);} //Bottom
		while(nodes[i].gety() > lengthD){nodes[i].sety(nodes[i].gety() - lengthD);} //Top
	
		if(nodes[i].getx() < 0 | nodes[i].getx() >lengthD | nodes[i].gety() < 0 | nodes[i].gety() > lengthD)
		{	
			cout << "ERROR 3" << endl;
			exit(1);
		}
	}
}

void Fibroblast::setFibroblastType(int fibroblastTypeD)
{
	fibroblastType = fibroblastTypeD;
}
int Fibroblast::getFibroblastType()
{
	return fibroblastType;
}		

double Fibroblast::getSpeed()
{
	return speed;
}
double Fibroblast::getRadiusC()
{
	return radiusC;
}

void Fibroblast::setAngleCI(double angleCID)
{
	angleCI = angleCID;
}
double Fibroblast::getAngleCI()
{
	return angleCI;
}
void Fibroblast::setNumber(int numberD)
{
	number = numberD;
}
int Fibroblast::getNumber()
{
	return number;
}
void Fibroblast::setNumberCI(int numberCID)
{
	numberCI = numberCID;
}
int Fibroblast::getNumberCI()
{
	return numberCI;
}
std::vector<Node> Fibroblast::getNodes()
{
	return nodes;
}

void Fibroblast::setBinaryContact(int contactD) // Can be useful as a test to check the fibroblasts were recognising contacts with each other
{
	binaryContact = contactD;
}
int Fibroblast::getBinaryContact()
{
	return binaryContact;
}

void Fibroblast::setMobility(float mobilityD)
{
	mobility = mobilityD;
}
float Fibroblast::getMobility()
{
	return mobility;
}
void Fibroblast::setOverUnder(std::string overUnderD)
{
	overUnder = overUnderD;
}
std::string Fibroblast::getOverUnder()
{
	//cout << "inside getOverUnder, overUnder = " << overUnder << endl;
	return overUnder;
}
void Fibroblast::setNucBarrier(int nucBarrierD)
{
	nucBarrier = nucBarrierD;
}
int Fibroblast::getNucBarrier()
{
	return nucBarrier;
}
void Fibroblast::setAttractionAmplitude(double aaD)
{
	attractionAmplitude = aaD;
}
double Fibroblast::getAttractionAmplitude()
{
	return attractionAmplitude;
}
void Fibroblast::setAttractionAngle(double angleD)
{
	attractionAngle = angleD;
}
double Fibroblast::getAttractionAngle()
{
	return attractionAngle;
}
//----------------------------------------------------------------------------
void Fibroblast::findCorrespondingBox(int fibreGridLengthD, int lengthD)
{	
	previousBox = box;
	double xCoord = nodes[0].getx();
	double yCoord = nodes[0].gety();
	
	double b = (double) fibreGridLengthD/lengthD; // determine size of gridBoxes

	int n0xGridBox = floor(b*xCoord);
	int n0yGridBox = floor(b*yCoord); 
	
	box = fibreGridLengthD*n0yGridBox + n0xGridBox; // working out which gridBox nH falls in 
}
void Fibroblast::setBox(int boxD)
{
	box = boxD;
}
int Fibroblast::getBox()
{
	return box;
}
int Fibroblast::getPreviousBox()
{
	return previousBox;
}

void Fibroblast::setTempNewAngle(double tempNewAngleD)
{
	tempNewAngle = tempNewAngleD;
}
double Fibroblast::getTempNewAngle()
{
	return tempNewAngle;
}
void Fibroblast::findVoxBox(int numberVoxBoxD, int lengthD)
{	
	double xCoord = nodes[0].getx();
	double yCoord = nodes[0].gety();
	
	double b = (double) numberVoxBoxD/lengthD; // determine size of gridBoxes
	int n0xGridBox = floor(b*xCoord);
	int n0yGridBox = floor(b*yCoord); 
	
	voxBox = numberVoxBoxD*n0yGridBox + n0xGridBox; // working out which gridBox nH falls in 	
}
int Fibroblast::getVoxBox()
{
	return voxBox;
}

void Fibroblast::setContactAngles(vector<double> contactAnglesD)
{
	contactAngles.clear();
	contactAngles = contactAnglesD;
}
void Fibroblast::setContactAnglesCC(vector<double> contactAnglesCCD)
{
	contactAnglesCC.clear();
	contactAnglesCC = contactAnglesCCD;
}
std::vector<double> Fibroblast::getContactAngles()
{
	return contactAngles;
}
std::vector<double> Fibroblast::getContactAnglesCC()
{
	return contactAnglesCC;
}
void Fibroblast::setCurrentContacts(std::vector<int> contactIDsD)
{
	currentContacts.clear();
	for(int i=0; i!=contactIDsD.size(); i++)
	{
		currentContacts.push_back(contactIDsD[i]);
	}

}
std::vector<int> Fibroblast::getCurrentContacts()
{
	return currentContacts;
}
void Fibroblast::setPreviousContacts()
{
	previousContacts.clear();
	for(int i=0; i!=currentContacts.size(); i++)
	{
		previousContacts.push_back(currentContacts[i]);
	}
}
std::vector<int> Fibroblast::getPreviousContacts()
{
	return previousContacts;
}
