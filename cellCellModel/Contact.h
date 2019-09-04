//9/5/16
#ifndef CONTACT_H
#define CONTACT_H
#include<string>
#include "Node.h"

class Contact
{	
	public:
		Contact();
		~Contact();	
		
		// when a contact happens we need to record all these things
			
		void setNodeTypeI(std::string nodeTypeD); // the node of cell i involved in the collision
		void setNodeTypeJ(std::string nodeTypeD);
		void setAngleI(double angleID); // the angle which cell i had been travelling at when it collided with cell j
		void setAngleJ(double angleJD);
		void setAngleDifference(); // the collision angle
		void setCellI(int CellD); // the identification number of cell i
		void setCellJ(int CellD);
		
		std::string getNodeTypeI();
		std::string getNodeTypeJ();
		double getAngleI();
		double getAngleJ();
		double getAngleDifference();
		int getCellI();
		int getCellJ();
				
	private:
		std::string nodeTypeI, nodeTypeJ; // H,S,B,T stands for Head,Shoulders, Body, Tail. 
		int cellI, cellJ; // identification number of the cells
		double angleI, angleJ, angleDifference; 
};

#endif

