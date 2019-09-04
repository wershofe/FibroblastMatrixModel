// 19/10/16
#ifndef CC_H
#define CC_H
#include "Contact.h"
#include<vector>
#include<string>
#include "Node.h"

class CC
{	
	public:
		CC(double radiusCD);
		~CC();		
	
		void setRadiusC(double radiusCD); // default mean radius of each node. 
		
		void setFirstPosition(int lengthD, double radiusCD, double thetaD); // random position and orientation
		void setFirstPosition(double xD, double yD, double radiusCD, double thetaD); // set position, random orientation - useful for making spheroid assays

		void setAngleCI(double angleCID); // angle of change triggered by CIL
		double getAngleCI(); // at beginning of each cycle this is set to 0, if CIL is triggered it gets changed

		double getRadiusC();
		
		void setNumber(int numberD); // give each cc a fixed number for identification
		int getNumber();
		void setNumberCI(int numberCID); // mark the number of the cc j which triggered the CI in cc i
		int getNumberCI();
		std::vector<Node> getNodes();

		void setBinaryContact(int contactD); // used as a test, to check ccs have identified contacts
		int getBinaryContact();
		double getTheta();
				
	private:
		double radiusC, theta; 
		bool mobility;
		std::vector<Node> nodes; // H,S,B,T stands for Head,Shoulders, Body, Tail. 
		std::vector<double> originalPosition; // useful if wanting to compute directionality ratio
		int binaryContact; // for testing purposes - has cc recognised collision
		double angleCI; // when CIL is activated, what is the angle response.
		int number; // each cc has a fixed number for identification
		int numberCI; // the number of the cc j that has activated the CI in cc i
};

#endif

