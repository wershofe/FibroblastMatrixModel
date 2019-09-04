// 19/10/16
#ifndef FIBROBLAST_H
#define FIBROBLAST_H
#include "Contact.h"
#include<vector>
#include<string>
#include "Node.h"

class Fibroblast
{	
	public:
		Fibroblast(double radiusCD, double speedMeanD, double speedSDD);
		~Fibroblast();		
		void setSpeed(double speedD); // fixed speed of fibroblast
		void setSpeed(double speedMeanD, double speedSDD); // for now speed follows a gaussian distribution and remains constant for each fibroblast
		void setRadiusC(double radiusCD); // default mean radius of each node. 
		
		void setFirstPosition(int lengthD, double radiusCD, int aspectRatioD); // random position and orientation
		void setFirstPosition(double xD, double yD, double radiusCD, int aspectRatioD); // set position, random orientation - useful for making spheroid assays
		void setFirstPosition(double xD, double yD, double thetaD, double radiusCD, int aspectRatioD); // set position and orientation - test specific collisions
		
		void setAngle(double angleD); // random angle
		void setPosition(double radiusCD, int lengthD, int aspectRatioD); // updates position of each fibroblast
		void applyBoundaryConditions(int lengthD); 

		void setAngleCI(double angleCID); // angle of change triggered by CIL
		double getAngleCI(); // at beginning of each cycle this is set to 0, if CIL is triggered it gets changed

		double getSpeed();
		double getTheta();
		double getRadiusC();
		
		void setNumber(int numberD); // give each fibroblast a fixed number for identification
		int getNumber();
		void setNumberCI(int numberCID); // mark the number of the fibroblast j which triggered the CI in fibroblast i
		int getNumberCI();
		std::vector<Node> getNodes();
		
		void setFibroblastType(int fibroblastTypeD);
		int getFibroblastType();
		
		void setBinaryContact(int contactD); // used as a test, to check fibroblasts have identified contacts
		int getBinaryContact();
		
		void setMobility(float mobilityD);
		float getMobility();
		
		void setOverUnder(std::string overUnderD);
		std::string getOverUnder();
		void setNucBarrier(int nucBarrierD);
		int getNucBarrier();
		
		void setAttractionAmplitude(double aaD);
		double getAttractionAmplitude();
		void setAttractionAngle(double angleD);
		double getAttractionAngle();
		
		void findCorrespondingBox(int fibreGridLengthD, int lengthD);
		void setBox(int boxD);
		int getBox();	
		int getPreviousBox();	
		
		void setTempNewAngle(double tempNewAngleD);
		double getTempNewAngle();
		
		void findVoxBox(int numberVoxBoxD, int lengthD);
		int getVoxBox();
		
		void setContactAngles(std::vector<double> contactAnglesD);
		void setContactAnglesCC(std::vector<double> contactAnglesCCD);
		std::vector<double> getContactAngles();
		std::vector<double> getContactAnglesCC();
		
		void setCurrentContacts(std::vector<int> contactIDsD);
		std::vector<int> getCurrentContacts();
		void setPreviousContacts();
		std::vector<int> getPreviousContacts();
	private:
		double speed, theta, radiusC; 
		float mobility;
		std::vector<Node> nodes; // H,S,B,T stands for Head,Shoulders, Body, Tail. 
		std::vector<double> originalPosition; // useful if wanting to compute directionality ratio
		int binaryContact; // for testing purposes - has fibroblast recognised collision
		double angleCI; // when CIL is activated, what is the angle response.
		int number; // each fibroblast has a fixed number for identification
		int numberCI; // the number of the fibroblast j that has activated the CI in fibroblast i
		int fibroblastType;
		std::string overUnder;
		int nucBarrier;
		int box; // for ECM
		int previousBox; // need to deposit fibres in previous box
		
		double attractionAmplitude; // 8-6-17 source chemoattractant 
		double attractionAngle;
		
		double tempNewAngle; // for parallelisation = computeAngle()...
		int voxBox; //for spatially dividing up
		
		std::vector<double> contactAngles;
		std::vector<double> contactAnglesCC;
		std::vector<int> previousContacts;
		std::vector<int> currentContacts;
};

#endif

