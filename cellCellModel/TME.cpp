// 19/10/16
#include <iostream>
#include "TME.h"
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <stdio.h>
#include "Fibroblast.h"
#include "Contact.h"
#include <vector>
#include <map>
#include <string>
#include <string.h>
#include <algorithm>
#include "Parameter.h"
#include "FibreBox.h"
#include <sstream>
#include <fstream>
#include <omp.h>

using namespace std;

TME::TME()
{
}

TME::~TME()
{
}

double TME::computeDistance2(Fibroblast cD, Fibroblast bD) // Euclidean distance but considering periodic boundary conditions
{
	double x1p = cD.getNodes()[0].getx(); //same coords with periodic boundaries applied
	double x2p = bD.getNodes()[0].getx();
	double y1p = cD.getNodes()[0].gety();
	double y2p = bD.getNodes()[0].gety();
	
	double dx, dy, dist;
	boundaryChecks(x1p,x2p,y1p,y2p);	

	dx = abs(x1p - x2p);
	dx = min(dx, lengthD-dx);
	dy = abs(y1p - y2p);
	dy = min(dy, lengthD-dy);
	
	if(dx<0 )
	{
		cout << "ERROR ASPARAGUS1" << endl;
		exit(1); // new
	}
	if(dy<0)
	{	
		cout << "ERROR ASPARAGUS2" << endl;
		exit(1); //new
	}
	if(dx>lengthD)
	{
		cout << "ERROR ASPARAGUS3" << endl;
		exit(1); //new
	}
	if(dy>lengthD)
	{
		cout << "ERROR ASPARAGUS4" << endl;
		exit(1); //new
	}
	dist = dx*dx + dy*dy;
	return dist;
}

void TME::boundaryChecks(double &x1D, double &x2D, double &y1D, double &y2D)
{
	while(x1D>lengthD)
	{
		x1D -= lengthD;
	}
	while(x2D>lengthD)
	{
		x2D -= lengthD;
	}
	while(y1D>lengthD)
	{
		y1D -= lengthD;
	}
	while(y2D>lengthD)
	{
		y2D -= lengthD;
	}
	
	while(x1D<0)
	{
		x1D += lengthD;
	}
	while(x2D<0)
	{
		x2D += lengthD;
	}
	while(y1D<0)
	{
		y1D += lengthD;
	}
	while(y2D<0)
	{
		y2D += lengthD;
	}
}

double TME::computeDistance2(double x1, double y1, double x2, double y2) // Euclidean distance but considering periodic boundary conditions
{
	double x1p = x1; //same coords with periodic boundaries applied
	double x2p = x2;
	double y1p = y1;
	double y2p = y2;
	
	double dx, dy, dist;
	boundaryChecks(x1p,x2p,y1p,y2p);
	
	dx = abs(x1p - x2p);
	dx = min(dx, lengthD-dx);
	dy = abs(y1p - y2p);
	dy = min(dy, lengthD-dy);
	
	if(dx<0 )
	{
		cout << "ERROR CHUTNEY1" << endl;
		exit(1); // new
	}
	if(dy<0)
	{
		cout << "ERROR CHUTNEY2" << endl;
		exit(1); // new
	}
	if(dx>lengthD)
	{
		cout << "ERROR CHUTNEY3" << endl;
		exit(1); // new
	}
	if(dy>lengthD)
	{
		cout << "ERROR CHUTNEY4" << endl;
		exit(1); // new
	}
	
	dist = dx*dx + dy*dy;
	return dist;
}

void TME::initialiseTME()
{
	initialiseCCs();	// CCs come first because fibroblasts initialise around CCs
	initialiseFibreBoxes();	
	initialiseFibroblasts();
}

void TME::initialiseFibroblasts()
{
	frame = 0;
	vector<Fibroblast> fibroblastList;
	for(int i = 0; i < rhoD; i++)
	{
		Fibroblast c(radiusCD, speedMeanD, speedSDD);
		/* // setting specific starting conditions to test code
		if(i==0)
		{
			c.setFirstPosition(500, 500 ,  0, radiusCD, aspectRatioD);
		}
		if(i==1)
		{
			c.setFirstPosition(500+70, 500 + 50,  M_PI/2, radiusCD, aspectRatioD);
			c.setSpeed(0);
			
		}
		if(i==2)
		{
			c.setFirstPosition(500+20, 500 + 5 ,  0, radiusCD, aspectRatioD);	
		}*/

		//--------------------------------------------------------------------------------------------------
		// putting fibroblasts around tumor - should comment out if there are no CCs
		if(rhoCCD!=0)
		{
			bool overlap = true;
			int overlapNum=0;
			while(overlap == true)
			{ 
				overlapNum =0;
				c.setFirstPosition(lengthD, radiusCD, aspectRatioD); // all random positions	
				double x1 = c.getNodes()[0].getx();
				double y1 = c.getNodes()[0].gety();
				for(int j=0; j!=CCs.size(); j++)
				{
					double x2 = CCs[j].getNodes()[0].getx();
					double y2 = CCs[j].getNodes()[0].gety();
					if(computeDistance2(x1,y1,x2,y2)<((radiusCCD+radiusCD)*(radiusCCD+radiusCD))) // if overlapping with CCs then no good
					{
						overlapNum++;
						break;
					}
					// if tumour is ring-shaped and want to initalise CAFs outside of tumour
				}
				if((x1-sourceXD)*(x1-sourceXD) + (y1-sourceYD)*(y1-sourceYD) <20000) // if within 100000 of CCs then no good
				{
					overlapNum++; 
					cout << "ERROR LAMB" << endl;
					exit(-1);
				} 
				if(overlapNum==0)
				{
					overlap=false;
				}
			}
		}
		//---------------------------------------------------------------------------------------------------
		if(woundSizeD!=0)
		{
			bool overlap = true;
			while(overlap == true)
			{ 
				c.setFirstPosition(lengthD, radiusCD, aspectRatioD); // all random positions	
				double x1 = c.getNodes()[0].getx();
				double y1 = c.getNodes()[0].gety();
				
				if((x1-512)*(x1-512) + (y1-512)*(y1-512) > woundSizeD*woundSizeD) // if within 100000 of CCs then no good
				{
					overlap=false;
				}
			}
		}
		
		//-----------------------------------------------------------------------------------------------------------
		else
		{
			c.setFirstPosition(lengthD, radiusCD, aspectRatioD); // all random positions
		}
		//c.setAngle((10*i)%30); // initial distribution quite aligned
		//c.setFirstPosition(400+(10*i),400+(10*i), 0.1*i, radiusCD); // specific position
		c.setNumber(fibroblastList.size()); // since we haven't yet pushed this fibroblast onto the fibroblasts vector
		c.findCorrespondingBox(fibreGridLengthD, lengthD); // set to find corresponding ECM box
		c.findVoxBox(numberVoxBoxD, lengthD); // set to find corresponding ECM voxBox
		fibroblastList.push_back(c);
	}
	
	fibroblasts = fibroblastList;
	setBoxToFibs(); // spatial segmentation - find voxBox
	setBoxToVoxFibs(); // spatial segmentation - find adjacent voxBoxes
}

void TME::setBoxToFibs() // assigns each fibroblasts to a voxBox
{
	boxToFibs.clear();
	for(int i=0; i!=fibroblasts.size(); i++)
	{
		int box = fibroblasts[i].getVoxBox();
		boxToFibs[box].push_back(fibroblasts[i]);
	}
}
void TME::setBoxToVoxFibs() // rewrites map between voxBox number and all the fibroblasts in that box and adjacent boxes
{
	boxToVoxFibs.clear(); 
	for(int i=0; i!=numberVoxBoxD*numberVoxBoxD; i++)
	{	
		computeAdjacentBoxes(i); // work out all adjacent boxes to box i
	}	
}

void TME::addFibroblast(int &timeD)
{
	
	int dividing = rhoD * pow (2, timeD/proliferationD);
	while(proliferationD && (fibroblasts.size()<dividing) && (fibroblasts.size() < 3000)) // currently when this number of fibroblasts, stops proliferating - ie maximum density
	{ 
		vector<int> freeBoxes;
		for(int i=0; i!= numberVoxBoxD*numberVoxBoxD; i++)
		{
			if(boxToFibs[i].size()==0)
			{
				freeBoxes.push_back(i);
			}
		}
		if(freeBoxes.size() !=0)
		{
			random_shuffle(freeBoxes.begin(),freeBoxes.end());
			for(int i=0; i!= freeBoxes.size(); i++)
			{
				int k=freeBoxes[i];
				if(boxToVoxFibs[k].size()!=0)
				{
					vector<int> pickFib;
					for(int j=0; j!=boxToVoxFibs[k].size(); j++)
					{
						pickFib.push_back(j);
					}
					random_shuffle(pickFib.begin(),pickFib.end());
					double theta = boxToVoxFibs[k][pickFib[0]].getTheta();
					int y = floor(k/numberVoxBoxD);
					int x = k - (numberVoxBoxD*y);
					double sizeOfFreeBox = lengthD/numberVoxBoxD;
					double rand01 = ((double) rand() / (RAND_MAX))*sizeOfFreeBox; //place new fibroblast randomly in that box
					double rand02 = ((double) rand() / (RAND_MAX))*sizeOfFreeBox;
					double yCoord = (y * sizeOfFreeBox) + rand01;
					double xCoord = (x * sizeOfFreeBox) + rand02;	
		
					Fibroblast c(radiusCD, speedMeanD, speedSDD);
					c.setFirstPosition(xCoord, yCoord, theta+M_PI, radiusCD, aspectRatioD);
					c.setNumber(fibroblasts.size()); // since we haven't yet pushed this fibroblast onto the fibroblasts vector
					c.findCorrespondingBox(fibreGridLengthD, lengthD); // set to find corresponding ECM box
					c.findVoxBox(numberVoxBoxD, lengthD); // set to find corresponding ECM voxBox
					fibroblasts.push_back(c);
					break;
				}
			}
		}
	}
}
void TME::computeAdjacentBoxes(int boxD) // want to know all fibroblasts in a voxel box and adjacent boxes - only those fibroblasts can be touching
{
	int x = boxD% numberVoxBoxD; // work out what the x and y markers are for that box index
	int y = (boxD-x)/numberVoxBoxD;
	
	int b = y-1, t = y+1, l = x-1, r = x+1; // working out adjacent boxes
    		if(b<0){b+=numberVoxBoxD;}
    		if(l<0){l+=numberVoxBoxD;}
    		if(t>=numberVoxBoxD){t-=numberVoxBoxD;} 
	   	if(r>=numberVoxBoxD){r-=numberVoxBoxD;}
	   	
	   	vector<int> adjacentFreeBox; // working out adjacent boxes
	   	int TL =  numberVoxBoxD*t + l; 
	   	int TM =  numberVoxBoxD*t + x; 
	   	int TR =  numberVoxBoxD*t + r; 
	   	int ML =  numberVoxBoxD*y + l; 
	   	int MR =  numberVoxBoxD*y + r; 
	   	int BL =  numberVoxBoxD*b + l; 
	   	int BM =  numberVoxBoxD*b + x; 
	   	int BR =  numberVoxBoxD*b + r; 

	   	adjacentFreeBox.push_back(TL);
	   	adjacentFreeBox.push_back(TM);
	   	adjacentFreeBox.push_back(TR);
	   	adjacentFreeBox.push_back(ML);
	   	adjacentFreeBox.push_back(boxD); // need to push back all cells in the same  box too!
	   	adjacentFreeBox.push_back(MR);
	   	adjacentFreeBox.push_back(BL);
	   	adjacentFreeBox.push_back(BM);
	   	adjacentFreeBox.push_back(BR);
	   
	   	for(int i=0; i!=adjacentFreeBox.size(); i++) // periodic boundary stuff with the boxes
		{
			while(adjacentFreeBox[i]>=numberVoxBoxD*numberVoxBoxD)
			{
				adjacentFreeBox[i]-= (numberVoxBoxD*numberVoxBoxD);
			}
		}
		
		for(int i=0; i!=adjacentFreeBox.size(); i++)
		{
			int adBox = adjacentFreeBox[i];
			for(int j=0; j!= boxToFibs[adBox].size(); j++) // vector with fibroblasts in all adjacent boxes (and same box)
			{
				boxToVoxFibs[boxD].push_back(boxToFibs[adBox][j]);
			}
		}
}


void TME::initialiseFibreBoxes() // if wanting to start with matrix defined in certain way
{
	if(IMTD =="A" || IMTD =="N")
	{
		for(int y=0; y!=fibreGridLengthD; y++)
		{
			for(int x=0; x!=fibreGridLengthD; x++)
			{
				int box = y*fibreGridLengthD + x;
			
				/*if(y==3) // eg if want to create a line of matrix across the top
				{
					FibreBox b(x,y,box,0,0,depRateD,reRateD,degRateD); 
					fibreBoxes.push_back(b);
				}
				else
				{
					FibreBox b(x,y,box,0,0,depRateD,reRateD,degRateD); // 
					fibreBoxes.push_back(b); // each fibreBox contains a record of the deposition, rearrangement and degradation rates
				}*/
				bool wounded = false;
				if(IMTD == "A")
				{
					FibreBox b(x,y,box,0,IMID,depRateD,reRateD,degRateD, noBinsD, wounded); 
					fibreBoxes.push_back(b); // initial aligning puts fibres in bin 0 always
				}
				else if(IMTD == "N")
				{
					double u1 = (double) rand() / (RAND_MAX);
					int bin = u1*noBinsD; // have 8 fibre bins per box
					FibreBox b(x,y,box,bin,IMID,depRateD,reRateD,degRateD, noBinsD, wounded); 
					fibreBoxes.push_back(b);
				}
				else
				{
					cout << "error HOT DOG" << endl;
					exit(35);
				}
			}
		}
	}
	else
	{
		string fileNameD = IMTD;
		// readInFile and assign boxes
		
		int box, x, y, density1, density2, density3, density4, density5, density6, density7, density8;	
		ifstream trajectoryFile(fileNameD.c_str());	
		
		while(trajectoryFile  >> box >> x >> y >> density1 >> density2 >> density3 >> density4 
								>> density5 >> density6 >> density7 >> density8 )
		{			
			vector<int> densities;
			densities.push_back(density1);
			densities.push_back(density2);
			densities.push_back(density3);
			densities.push_back(density4);
			densities.push_back(density5);
			densities.push_back(density6);
			densities.push_back(density7);
			densities.push_back(density8);
			
			bool wounded=woundMatrix(x,y);
			if(wounded==true)
			{
				densities.clear();
				for(int i=0; i!=8; i++)
				{
					densities.push_back(0);
				}
			}
			
			FibreBox b(x,y,box,densities,depRateD, reRateD, degRateD, noBinsD, wounded);
			fibreBoxes.push_back(b);
		}	
	}
}

bool TME::woundMatrix(int &x, int &y)
{
	bool wounded=false;
	
	double b = (double) fibreGridLengthD/lengthD; 
	double xCoord = x/b;
	double yCoord = y/b;
	double dist = computeDistance2(xCoord, yCoord, 512, 512);
	if(dist<woundSizeD*woundSizeD)
	{
		wounded=true;
	}
	
	return wounded;
}

void TME::initialiseCCs()
{
	double divideArc = M_PI/(rhoCCD*rhoCCD);
	double divideCircle = 2*M_PI/(rhoCCD*rhoCCD);
	for(int i=0; i!=rhoCCD; i++)
	{
		for(int j=0; j!=rhoCCD; j++)
		{
			CC c(radiusCCD);
			double segment = j+rhoCCD*i;
			
			// circle
			c.setFirstPosition(sourceXD+(10*radiusCCD*cos(divideCircle*segment)),sourceYD+(10*radiusCCD*sin(divideCircle*segment)),radiusCCD,
							divideCircle*segment + M_PI/2); // need to play with these numbers to get good shape
			
			// bay
			/*c.setFirstPosition(512+(8*radiusCCD*cos(divideArc*segment)),512+(8*radiusCCD*sin(divideArc*segment)),radiusCCD,
							divideArc*segment + M_PI/2);*/
			
			//rectangle
			/*c.setFirstPosition(512+(2*radiusCCD*i),512+(2*radiusCCD*j),radiusCCD);*/		
			
			c.setNumber(j + rhoCCD*i); // since we haven't yet pushed this fibroblast onto the fibroblasts vector
			CCs.push_back(c);
		}		
	}	
	//sink - if want to have a morphological instability
	/*
	CC c1(radiusCCD);
	c1.setFirstPosition(sourceXD,sourceYD+(11*radiusCCD),radiusCCD, M_PI/2); // need to play with this number to get good shape
	c1.setNumber(CCs.size());
	CCs.push_back(c1);
	
	CC c2(radiusCCD);
	c2.setFirstPosition(sourceXD,sourceYD+(12*radiusCCD),radiusCCD, M_PI/2); // need to play with this number to get good shape
	c2.setNumber(CCs.size());
	CCs.push_back(c2);
	*/
}

vector<Contact> TME::findContacts(Fibroblast &cD, int &timeD) // identifies potential contacts
{ // Find potential contacts, then will check them properly in checkContacts, then will put them into contacts vector
	vector<Contact> contacts;
	cD.setBinaryContact(0);
	vector<Fibroblast> adjFibs = boxToVoxFibs[cD.getVoxBox()]; // only consider potential contacts out of cells which are in adjacent voxBoxes
	
	for(int j=0; j!= adjFibs.size(); j++)
	{
		if(adjFibs[j].getNumber() != cD.getNumberCI() && adjFibs[j].getNumber() != cD.getNumber() ) 
		// makes fibroblast that triggered the CI invisible (analogous to repolarisation)
		{//only consider fibroblast j if it the heads are less than 14r (distance if they were tail-to-tail) apart (otherwise for sure not touching)
			double maxDist = computeDistance2(cD,adjFibs[j]);
			if(maxDist < ((aspectRatioD+2)*2*radiusCD)*((aspectRatioD+2)*2*radiusCD)) // check if it's possible for cells to be touching
			{	// look at all nodes of fibroblast i and fibroblast j
				
				for(int a=0; a < cD.getNodes().size(); a++) 
				{
					for(int b=0; b < adjFibs[j].getNodes().size(); b++)
					{
						checkContacts(cD, adjFibs[j], a, b, timeD, contacts);
					}
				}
			}	
		}
	}
	vector<int> contactIDs;
	for(int i=0; i!=contacts.size(); i++)
	{
		int j=contacts[i].getCellJ();
		contactIDs.push_back(j);
	}
	cD.setCurrentContacts(contactIDs);

	return contacts;
}

vector<Contact> TME::findContactsCC(Fibroblast &cD, int &timeD)
{ // Find potential contacts, then will check them properly in checkContacts, then will put them into contacts vector
	vector<Contact> contactsCC;
	cD.setBinaryContact(0);

	for(int j = 0; j!= CCs.size(); j++)
	{
		for(int a=0; a < cD.getNodes().size(); a++) 
		{
			checkContactsCC(cD, CCs[j], a, timeD, contactsCC);

		}
	}	
	return contactsCC;
}

void TME::checkContacts(Fibroblast &cD, Fibroblast &bD, int iD, int jD, int &timeD, vector<Contact> &contactsD) // identifies true contacts
{		
	string c, b;
	c = cD.getNodes()[iD].getNodeType();
	b = bD.getNodes()[jD].getNodeType();
	
	if(c!="C" && b!="C") // C doesn't exist as a node
	{
		double touching; // works out the distance between body parts for a touch to occur - it will be smaller if a head is involved
		if( c == "H" && b == "H")
		{
			//touching = radiusCD; // reminder: radiusCD is the radius of the larger (body) nodes
			touching = (2*headSizeD)*radiusCD;
		}
		else if( (c == "H" && b == "B") ||  (c == "B" && b == "H") )
		{
			touching = (headSizeD+1)*radiusCD;
		}
		else if( (c == "H" && b == "T") ||  (c == "T" && b == "H") )
		{
			touching = (headSizeD+0.5)*radiusCD;
		}
		else if( c == "B" && b == "B" )
		{
			touching = (1+1)*radiusCD;
		}
		else if( (c == "B" && b == "T") ||  (c == "T" && b == "B") )
		{
			touching = (1+0.5)*radiusCD;
		}
		else if( c == "T" && b == "T" )
		{
			touching = (0.5+0.5)*radiusCD;
		}
		
		else
		{
			cout << "error pancake c b = " <<c << "\t" << b << endl;
			exit(24);
		}
			double distance = computeDistance2(cD.getNodes()[iD].getx(), cD.getNodes()[iD].gety(), 
						bD.getNodes()[jD].getx(), bD.getNodes()[jD].gety());	
						
			//cout <<  cD.getNodes()[iD].getNodeType() << ": (" << cD.getNodes()[iD].getx() << "\t" << cD.getNodes()[iD].gety() << "), " << 
			//bD.getNodes()[jD].getNodeType() << ": (" << 
			//bD.getNodes()[jD].getx() << "\t" << bD.getNodes()[jD].gety() << ") "
			//<<": distance touching*touching = " << distance << "\t" << touching*touching << endl;
			if( distance <= touching*touching) 
			{	// if fibroblasts are touching, a contact is registered and added to the vector of contacts for fibroblast i
			//	cout << "CONTACT!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
				Contact con;
				con.setNodeTypeI(cD.getNodes()[iD].getNodeType()); // indicate the node of fibroblast i that was involved in the collision

				con.setNodeTypeJ(bD.getNodes()[jD].getNodeType());
				con.setCellI(cD.getNumber());
				con.setCellJ(bD.getNumber());
				con.setAngleI(cD.getTheta());
				con.setAngleJ(bD.getTheta());
				con.setAngleDifference();

				contactsD.push_back(con);
				if(CILaD<500)
				{
					checkCIL(cD, bD, iD, jD, con, timeD); // need to check if that contact is triggering CIL or not
				}
			
			}
	}
}

void TME::freezeCell(int numCell) // make the cell freeze
{
	if(fibroblasts[numCell].getNumber() != numCell)
	{
		cout << "fibroblasts[numCell].getNumber(), numCell = " << fibroblasts[numCell].getNumber() << "\t" << numCell << endl;
		exit(-1);
	}
	fibroblasts[numCell].setMobility(0.0); 
	
}

void TME::checkContactsCC(Fibroblast &cD, CC  &ccD, int iD, int &timeD, vector<Contact> &contactsCCD)
{
	double touching; // works out the distance between body parts for a touch to occur - it will be smaller if a head is involved

	if(cD.getNodes()[iD].getNodeType()!="C")
	{
		if(cD.getNodes()[iD].getNodeType() == "H")
		{
			touching = 0.5*radiusCD + radiusCCD;
		}
		else
		{
			touching = radiusCD + radiusCCD;
		}
		
		if( computeDistance2(cD.getNodes()[iD].getx(), cD.getNodes()[iD].gety(), 
				ccD.getNodes()[0].getx(), ccD.getNodes()[0].gety()) < touching*touching) 
		{	// if fibroblasts are touching, a contact is registered and added to the vector of contacts for fibroblast i
			Contact con;
			con.setNodeTypeI(cD.getNodes()[iD].getNodeType()); // indicate node of fibroblast i that was involved in the collision

			con.setNodeTypeJ("CC");
			con.setCellI(cD.getNumber());
			con.setCellJ(ccD.getNumber());
			con.setAngleI(cD.getTheta());
			con.setAngleJ(ccD.getTheta());
			con.setAngleDifference();

			contactsCCD.push_back(con);
			
			checkHetCIL(cD, ccD, iD, con, timeD); // need to check if that contact is triggering CIL or not
		}
	}	
}

// Circumstances under which CIL is activated
void TME::checkCIL(Fibroblast &cD, Fibroblast &bD, int iD, int jD, Contact conD, int &timeD)
{
	bool CIL = false;
	bool alreadyTouching=false; // only want to consider doing CIL for a new contact
	for(int i=0; i!=cD.getPreviousContacts().size(); i++)
	{
		if(bD.getNumber()==cD.getPreviousContacts()[i])
		{
			alreadyTouching = true;
			break;
		}
	}
	
	double randNo = ((double) rand() / (RAND_MAX));
	if(alreadyTouching == false && randNo<probCILD)
	{
		double lowerBound, upperBound; // for size of pizza slices
		if(CILaD>=90) // basically when we want to turn off CIL if set CILa >= 90 (apart from in head-head collisions)
		{
			lowerBound = M_PI;
			upperBound = M_PI;
		}
		else
		{
			lowerBound = CILaD*(M_PI/180);
			upperBound = (180-CILaD)*(M_PI/180); // nematic alignment
		}	
		// checking whether or not to activate CIL depending on the nodes involved in the collision
		if(cD.getNodes()[iD].getNodeType() == "H" && bD.getNodes()[jD].getNodeType() == "H" ) // coded so CIL can only be triggered on HH collision
		{	
			if(conD.getAngleDifference() > CILaD * (M_PI/180))// 
			{
				CIL=true;
			}
		}
		// NB CIL cannot be triggered by collision into the tail (T) of fibroblast j
		else if(cD.getNodes()[iD].getNodeType() == "H" && (bD.getNodes()[jD].getNodeType() == "B") )
		{
			if(conD.getAngleDifference() > lowerBound && conD.getAngleDifference() < upperBound) // ie if CIL is triggered 
			{			
				CIL = true;
			}
		}
	}
	if(CIL==true)
	{
		double cilResponse = drawFromRose(); 
		
		double thetaJprime = bD.getTheta() - cD.getTheta();
		double gradient = sin(thetaJprime)/cos(thetaJprime);
		if(gradient<0)
		{
			cilResponse *= -1; // so that the cell i bounces of cell j the correct way to make it less criss-crossed
		} 
		cD.setAngleCI(cilResponse);
		cD.setNumberCI(bD.getNumber());
	}	
}

double TME::drawFromRose() // fibroblast's response to CIL can be fixed, or drawn from a distribution - 
											// so far the options are a normal or uniform distribution
{
	double adjustAngle = 0;
	
	if(CILtypeD == "A") // ~N(0,.1) - "A" stands for aligning
	{
		double u1 = (double) rand() / (RAND_MAX);
		double u2 = (double) rand() / (RAND_MAX);
		double z0 = fakeRoseMeanD+fakeRoseSDD*sqrt(-2*log(u1))*cos(2*M_PI*u2); // box-muller transform 
		adjustAngle = z0;
		
		while(adjustAngle>180)
		{
			adjustAngle-=360;
		}
		while(adjustAngle<-180)
		{
			adjustAngle+=360;
		}
		adjustAngle = adjustAngle*(M_PI/180);
	}
	else if(CILtypeD == "N") // ~U[0,2*M_PI] - "N" stands for non-aligning
	{
		adjustAngle = 2*M_PI*((double) rand() / (RAND_MAX)); 
	}
	else if(CILtypeD == "F") // F = fixed angle of change
	{
		adjustAngle = CILbD*(M_PI/180); // this is CILb in the main - ie the fixed response angle of fibroblast i when CIL is triggered
	}
	else
	{
		cout << "ERROR TOMATO drawFromRose" << endl;
		exit(1);
	}
	return adjustAngle;
}

double TME::computeAngle(Fibroblast &cD, vector<Contact> &contactsD, vector<Contact> &contactsCCD, int &timeD)
{ // this is the assembling of the model equation, determining the orientation of each fibroblast at each time step.
	double w1Component=0, w5Component=0;
	int noConsideredFibres = 0;
	computeW1Component(cD, contactsD, w1Component, timeD); // cell persistence
	if(w2D!=0)
	{
		computeW2Component(cD, contactsD, timeD); // contact guidance
	}
	//computeW3Component(cD, timeD); // How much orientation of fibroblast i is influenced by chemoattractant (CAT)
	if(w4D!=0)
	{
		computeW4Component(cD, contactsCCD, timeD); // influence of CCs
	}
	if(w5D!=0)
	{
		computeW5Component(cD, timeD, noConsideredFibres, w5Component); //ECM influence
	}
	double angle = determineNewAngle(cD, w1Component, w5Component, noConsideredFibres); // equation of model
	return angle;
}

void TME::computeW1Component(Fibroblast &cD, vector<Contact> &contactsD, double &w1ComponentD, int &timeD)
{	
	//cout << "inside computeW1Component, contactsD.size() cD.getAngleCI() = " << contactsD.size() << "\t" << cD.getAngleCI() << endl;
	double u1 = (double) rand() / (RAND_MAX);
	double u2 = (double) rand() / (RAND_MAX);
	double persistencePolarity = polarityD*(sqrt(-2*log(u1))*cos(2*M_PI*u2)); // box-muller transform
	
	vector<int> uniqueCellContacts;
	for(int i=0; i!=contactsD.size(); i++)
	{
		int j=contactsD[i].getCellJ();
		if(find(uniqueCellContacts.begin(), uniqueCellContacts.end(), j) != uniqueCellContacts.end())
		{	
		}
		else
		{
			uniqueCellContacts.push_back(j);
		}
	}
	/*if(contactsD.size()!=0)
	{
		cout << "Omsode cp,[itew1 cpm[pentn, contactsD.size() = " << contactsD.size() << endl;
		cout << "Cell = " << cD.getNumber() << endl;
		cout << "nodeI nodeJ " << endl;
		for(int i=0; i!=contactsD.size(); i++)
		{
			cout << contactsD[i].getNodeTypeI() << "\t" << contactsD[i].getNodeTypeJ() << endl;
		}
	}*/
	if(uniqueCellContacts.size() == 1) // only triggering CIL when cell has single contact
	{
		w1ComponentD = cD.getTheta() + cD.getAngleCI() + persistencePolarity ; // reminder: cD.getAngleCI()=0 <=> CIL was not triggered
	}
	else
	{
		w1ComponentD = cD.getTheta() + persistencePolarity ; // 
	}
	
}

void TME::computeW2Component(Fibroblast &cD, vector<Contact> &contactsD, int &timeD)
{
	double upperBound=M_PI, lowerBound=M_PI; //default values - ie if CILa>90, so that upperBound is always >= lowerBound
	if(CILaD<90)
	{
		lowerBound = CILaD*(M_PI/180);
		upperBound = (180-CILaD)*(M_PI/180);
	}
	vector<double> w2ComponentVec;
	for(int i=0; i!=contactsD.size(); i++)
	{
		double w2Element=0;
		checkHeadHeadcontacts(contactsD[i], w2Element, cD);
		checkOtherContacts(contactsD[i], w2Element, cD);
		modifyW2Element(contactsD[i],w2Element);
		w2ComponentVec.push_back(w2Element);
	}
	cD.setContactAngles(w2ComponentVec); 
}

void TME::modifyW2Element(Contact &contactD, double &w2ElementD) // NEW 13-2-17 make it so that w2 is never more than pi/2 away from theta_i
{
	while((contactD.getAngleI()-w2ElementD)*(contactD.getAngleI()-w2ElementD)>(M_PI/2)*(M_PI/2))
	{
		double incUp = w2ElementD + M_PI; 
		double incDown = w2ElementD - M_PI;
		if((contactD.getAngleI()-incUp)*(contactD.getAngleI()-incUp) < (contactD.getAngleI()-w2ElementD)*(contactD.getAngleI()-w2ElementD))
		{
			w2ElementD = incUp;
		}
		else
		{
			w2ElementD = incDown;
		}
	}
}

void TME::checkHeadHeadcontacts(Contact &contactD, double &w2ElementD, Fibroblast &cD)
{
	if(contactD.getCellJ()!=cD.getNumberCI())
	{
		if(contactD.getNodeTypeI()=="H" && contactD.getNodeTypeJ()=="H")
		{
			if(contactD.getAngleDifference()<=M_PI/2)
			{	
				w2ElementD += contactD.getAngleJ();
			}
			else // nematic behaviour
			{
				reverseAngle(contactD, w2ElementD);
			}
			
		}
	}
}
void TME::checkOtherContacts(Contact &contactD, double &w2ElementD, Fibroblast &cD)
{
	if(contactD.getCellJ()!=cD.getNumberCI())
	{
		if(!(contactD.getNodeTypeI()=="H" && contactD.getNodeTypeJ()=="H"))
		{
			if(contactD.getAngleDifference()<=M_PI/2)
			{	
				w2ElementD += contactD.getAngleJ();
			}
			else // nematic behaviour
			{
				reverseAngle(contactD, w2ElementD);
			}
			
		}
	}
}
//------------------------------------------------------------------------------------------------------------------------
void TME::computeW4Component(Fibroblast &cD, vector<Contact> &contactsCCD, int &timeD) // CAF-CC interactions
{
	double upperBound=M_PI, lowerBound=M_PI; //default values - ie if CILa>90, so that upperBound is always >= lowerBound
	if(CILaD<90)
	{
		lowerBound = CILaD*(M_PI/180);
		upperBound = (180-CILaD)*(M_PI/180);
	}
	vector<double> w4ComponentVec;
	for(int i=0; i!=contactsCCD.size(); i++)
	{
		double w4Element=0;
		checkHetNematic(contactsCCD[i], w4Element, lowerBound, upperBound);
		modifyW4Element(contactsCCD[i],w4Element);
		w4ComponentVec.push_back(w4Element);
	}
	cD.setContactAnglesCC(w4ComponentVec);
}

void TME::modifyW4Element(Contact &contactCCD, double &w4ElementD) // NEW 13-2-17 make it so that w4 is never more than pi/2 away from theta_i
{
	while((contactCCD.getAngleI()-w4ElementD)*(contactCCD.getAngleI()-w4ElementD)>(M_PI/2)*(M_PI/2))
	{
		double incUp = w4ElementD + M_PI; 
		double incDown = w4ElementD - M_PI;
		if((contactCCD.getAngleI()-incUp)*(contactCCD.getAngleI()-incUp) < (contactCCD.getAngleI()-w4ElementD)*(contactCCD.getAngleI()-w4ElementD))
		{
			w4ElementD = incUp;
		}
		else
		{
			w4ElementD = incDown;
		}
	}
}

void TME::checkHetCIL(Fibroblast &cD, CC &bD, int iD, Contact conD, int &timeD)
{
	double cilResponse = sampleCIL(); // the angle of change triggered by CIL. This can be fixed or drawn from a distribution
	// this bit works out whether the fibroblast should rotate clockwise or acw, in order to be less "criss-crossing" with the passive fibroblast.

	double add = cD.getTheta() + cilResponse;
	double sub = cD.getTheta() - cilResponse;
	
	double a = abs(add - bD.getTheta());
	double a1 = min(a,(2*M_PI)-a);
	double a2 = min(a1, M_PI-a1);
	
	double s = abs(sub - bD.getTheta());
	double s1 = min(s,(2*M_PI)-s);
	double s2 = min(s1, M_PI-s1);
	
	if(s2<a2)
	{
		cilResponse = -cilResponse;
	}
	
	double lowerBound = M_PI, upperBound = M_PI; // for size of pizza slices
	if(hetCILaD<90)
	{
		lowerBound = hetCILaD*(M_PI/180);
		upperBound = (180-hetCILaD)*(M_PI/180);
	}	
	// checking whether or not to activate CIL depending on the nodes involved in the collision - for now coded so only potential for CIL when head collision
	if(cD.getNodes()[iD].getNodeType() == "H")
	{
		if(conD.getAngleDifference() > hetCILaD * (M_PI/180))// 
		{
			cD.setAngleCI(cilResponse);
			cD.setNumberCI(bD.getNumber());
		}
	}
}
void TME::checkHetNematic(Contact &contactCCD, double &w4ElementD, double &lowerBoundD, double &upperBoundD) // when CAFs hit CCs, nematic
{
	if(hetCILaD>180) // ie CIL is never triggered
	{
		if(contactCCD.getAngleDifference()<=M_PI/2)
		{	
			w4ElementD += contactCCD.getAngleJ();
		}
		else // nematic behaviour
		{
			reverseAngle(contactCCD, w4ElementD);
		}
	}
	else // head-head collisions only have left-hand of pizza wedge
	{
		if(contactCCD.getAngleDifference() <= hetCILaD*(M_PI/180) ) 
		{ 
			w4ElementD += contactCCD.getAngleJ();
		}
	}	
}
//------------------------------------------------------------------------------------------------------------------------
void TME::computeW5Component(Fibroblast &cD, int &timeD, int &noConsideredFibres, double &w5ComponentD) // CAF-ECM feedback
{	
	
	int box = cD.getBox(); 
	double theta = cD.getTheta(); 
	
	noConsideredFibres = fibreBoxes[box].getDominantDensity(); 
	double orientation = fibreBoxes[box].getDominantOrientation(); 
	double angleDifference = abs(theta-orientation);
	
	while(angleDifference > M_PI)
	{
		angleDifference = (2*M_PI) - angleDifference;
	}
	if(angleDifference > M_PI/2) // nematic - cells can move either way along a fibre
	{
		orientation += M_PI;
	}	
	
	w5ComponentD = orientation; 
}

//------------------------------------------------------------------------------------------------------------------------
void TME::reverseAngle(Contact &contactD, double &w2ElementD) // for nematic gliding behaviour
{
	double proposalAngle = contactD.getAngleJ() + M_PI;
	while(proposalAngle > M_PI) // NEW 12-2-17
	{
		proposalAngle -= 2*M_PI; // NEW 12-2-17
	} // NEW 12-2-17

	double reverseAngle = proposalAngle;
	w2ElementD += reverseAngle;
}

void TME::blockFromCCs() // stops CAFs going inside tumor
{
	for(int i=0; i!=fibroblasts.size(); i++)
	{
		for(int j=0; j!=CCs.size(); j++)
		{
			double x1 = fibroblasts[i].getNodes()[0].getx();
			double y1 = fibroblasts[i].getNodes()[0].gety();
			double x2 = CCs[j].getNodes()[0].getx();
			double y2 = CCs[j].getNodes()[0].gety();
			if(computeDistance2(x1,y1,x2,y2)<(0.25*(radiusCCD+radiusCD)*(radiusCCD+radiusCD)) ) //0.25 because allow for some contact guidance and movement
			{	
				fibroblasts[i].setMobility(0.0);
			}
		}
	}
}

double TME::computeLevelCAT(Fibroblast &cD)
{
	double x=cD.getNodes()[0].getx();
	double y=cD.getNodes()[0].gety();
	double b = (double)20/lengthD;
	
	int iBox = floor(b*x);
	int jBox = floor(b*y);
	
	int index = iBox + jBox*20;
	//cout << "b lengthD, iBox, jBox = " << b << "\t" << lengthD << "\t" << iBox << "\t" << jBox << endl;
	//cout << "cD.getNumber cD.getX index Mij[index] = " << cD.getNumber() << "\t" << x << "\t" << index << "\t" << Mij[index] << endl;
	return Mij[index];
}

void TME::computeW3Component(Fibroblast &cD, int &timeD)
{
	double distanceFromSource = computeDistance2(cD.getNodes()[1].getx(),cD.getNodes()[1].gety(),sourceXD,sourceYD);
	//cout << "Distance between (" << cD.getNodes()[1].getx() << "," << cD.getNodes()[1].gety() << ") and (" << (lengthD/2)
	//	<< "," << (lengthD/2) << ") is " << distanceFromSource << endl;
	double aa = computeLevelCAT(cD);
	cD.setAttractionAmplitude(aa); // for now a uniform area where attractant works
	
	// establish which quadrant
	double x1 = cD.getNodes()[1].getx();
	double y1 = cD.getNodes()[1].gety();
	double x2 = sourceXD;
	double y2 = sourceYD;
	
	double phi = atan((y2-y1)/(x2-x1));
	double pullAngle;
	if(x1>=x2 && y1>=y2) // quad 1
	{
		pullAngle = phi + M_PI;

	}
	if(x1<=x2 && y1>=y2) // quad 2
	{
		pullAngle = phi;

	}
	if(x1<=x2 && y1<=y2) // quad 3
	{
		pullAngle = phi;
		
	}
	if(x1>=x2 && y1<=y2) // quad 4
	{
		pullAngle = -(phi + M_PI);

	}	

	// establish whether cell should move cw or acw towards CAT
	
	double acwD = cD.getTheta() + (w3D*pullAngle);
	double cwD = cD.getTheta() - (w3D*pullAngle);
	cD.setAttractionAngle(pullAngle); // default wants to turn cell acw
	
	if( (sin(acwD)-sin(pullAngle))*(sin(acwD)-sin(pullAngle)) > (sin(cwD)-sin(pullAngle))*(sin(cwD)-sin(pullAngle)) )
	{
		cD.setAttractionAngle(-pullAngle); // when it is more efficient to turn cell cw
	}	
	else
	{
	}
}

double TME::determineNewAngle(Fibroblast &cD, double &w1ComponentD, double &w5ComponentD, int &noConsideredFibresD) // main equation of model
{	
	
	int diracW2 = 0; // if there are some cell contacts make nonzero
	int diracW5 = 0; // if there is some ECM below

	double X_w1 = cos(w1ComponentD);	
	double Y_w1 = sin(w1ComponentD);
	
	
	//-----------------------------------------------------
	
	double X_w2 = 0;
	double Y_w2 = 0;

	for(int i=0; i!=cD.getContactAngles().size(); i++)
	{
		X_w2 += cos(cD.getContactAngles()[i]);
		Y_w2 += sin(cD.getContactAngles()[i]);
	}
	if(cD.getContactAngles().size()!=0)
	{
		X_w2 /= cD.getContactAngles().size();
		Y_w2 /= cD.getContactAngles().size();
		diracW2 = 1;
	}
	
	//-----------------------------------------------------

	double X_w5 = cos(w5ComponentD);	
	double Y_w5 = sin(w5ComponentD);
	
	
	if(noConsideredFibresD<5)
	{
		diracW5 = 0;
	}
	else if(noConsideredFibresD>=5 && noConsideredFibresD<=10)
	{
		diracW5 = 0.2*(noConsideredFibresD-5);
	}
	else
	{
		diracW5 = 1;
	}
	
	
	
	/*
	if(noConsideredFibresD!=0)
	{
		diracW5 = 1;
	}
	*/
	//-----------------------------------------------------
	double X_componentNum = w1D*X_w1  + diracW2*w2D*X_w2 + diracW5*w5D*X_w5;
	double Y_componentNum = w1D*Y_w1  + diracW2*w2D*Y_w2 + diracW5*w5D*Y_w5;
	
	
	double denominator = w1D + diracW2*w2D + diracW5*w5D;
		
	double X_component = X_componentNum/denominator;
	double Y_component = Y_componentNum/denominator;		
	
	//-----------------------------------------------------	
	double dummyAngle = atan(Y_component/X_component);
	double angle;
	//arctan gives angle in range [-pi/2,pi/2] want to shift to range [0,2pi] depending on quadrant
	if(X_component>=0 && Y_component>=0) 
	{
		angle = dummyAngle;
	}
	if(X_component<0 && Y_component>=0)
	{
		angle = M_PI + dummyAngle;
	}
	if(X_component<0 && Y_component<0)
	{
		angle = M_PI + dummyAngle;
	}
	if(X_component>=0 && Y_component<0)
	{
		angle = 2*M_PI + dummyAngle;
	}
	return angle;
}

// This function is necessary in order to compute all the average angles before updatted any of the positions so updattes happen simultaneously.
vector<double> TME::assembleNewAngles(int &timeD) 
{	
	vector<double> newAngles;
	//omp_set_num_threads(4); // useful if running on farm
	// openMP parallelization!
	//# pragma omp parallel for 
	for(int i = 0; i < fibroblasts.size(); i++)
	{
		vector<Contact> contacts = findContacts(fibroblasts[i], timeD); 
		vector<Contact> contactsCC = findContactsCC(fibroblasts[i],timeD);
		double tempNewAngle = computeAngle(fibroblasts[i], contacts, contactsCC, timeD); // store tempNewAngle so that for loop can be parallelized
		fibroblasts[i].setTempNewAngle(tempNewAngle);
	}
	for(int i=0; i!=fibroblasts.size(); i++)
	{
		newAngles.push_back(fibroblasts[i].getTempNewAngle());
	}
	return newAngles;	
}

std::vector<CC> TME::getCCs()
{
	return CCs;
}

int TME::getNumberVoxBox()
{
	return numberVoxBoxD;
}
// the main updating function containing everything
void TME::updateTME(int &timeD)
{
	Mij.clear();
	readMatrix(timeD);
	proposals.clear();
	nucHeadContacts.clear();
	nucNucContacts.clear(); 
	for(int j=0; j!=fibroblasts.size(); j++)
	{
		fibroblasts[j].setPreviousContacts();
		fibroblasts[j].setMobility(1);
		fibroblasts[j].setAngleCI(0); // at the beginning of each run, reset CIL response of each fibroblast to 0
	}
	blockFromCCs();
	vector<double> angleList = assembleNewAngles(timeD); 
	proposals = fibroblasts;  
	for(int j = 0; j < proposals.size(); j++) // set new positions of the fibroblasts depending on their new orientations
	{
		//cout << "\t \t cell = " << j << " thread = " << id << endl;
		proposals[j].setAngle(angleList[j]);
		proposals[j].setPosition(radiusCD, lengthD, aspectRatioD);
		proposals[j].applyBoundaryConditions(lengthD);
	} 
	//stopOverlaps(); 
	if(speedMeanD!=0)
	{
		judderOverlaps();
	}
	for(int j = 0; j<fibroblasts.size(); j++) // set new positions of the fibroblasts depending on their new orientations
	{
		
		fibroblasts[j].setAngle(angleList[j]); 
		fibroblasts[j].setPosition(radiusCD, lengthD, aspectRatioD); 
		fibroblasts[j].applyBoundaryConditions(lengthD);
	} 
	frame++;
	timeD ++;
	/*
	if(frame%20==0)
	{
		//cout << frame << ", jamRatio = " << jamRatio << endl;
	}
	*/
	for(int j = 0; j != fibroblasts.size(); j++)
	{
		fibroblasts[j].findCorrespondingBox(fibreGridLengthD, lengthD);
		fibroblasts[j].findVoxBox(numberVoxBoxD, lengthD);
		depositFibres(fibroblasts[j]);
	}	
	setBoxToFibs();	
	setBoxToVoxFibs(); 
	addFibroblast(timeD);
}

void TME::judderOverlaps()
{
	for(int i=0; i<proposals.size(); i++)
	{	
		vector<Fibroblast> adjFibs = boxToVoxFibs[proposals[i].getVoxBox()]; // only checks volume exclusion for cells nearby
		
		for(int j=0; j!=adjFibs.size(); j++)
		{ 
			if(adjFibs[j].getNumber() != proposals[i].getNumber()) 
			{
				double maxDist = computeDistance2(proposals[i], adjFibs[j]);
				if(maxDist < ((aspectRatioD+2)*2*radiusCD)*((aspectRatioD+2)*2*radiusCD))
				{	
					checkVE(proposals[i], adjFibs[j]);
					
					//proposals[adjFibs[j].getNumber()].setOverUnder(adjFibs[j].getOverUnder());
					//proposals[adjFibs[j].getNumber()].setNucBarrier(adjFibs[j].getNucBarrier());						
				}	
			}
		}
	}
}

void TME::checkVE(Fibroblast &cD, Fibroblast &bD)
{	
	double x1,y1,x2,y2;
	
	int cellID = cD.getNumber();
	if(proposals[cellID].getNumber() != cellID || fibroblasts[cellID].getNumber() != cellID)
	{
		cout << "issue with cellID" << endl;
		exit(-1);
	}	
	
	for(int i=0; i!=bD.getNodes().size(); i++)
	{
		x1 = cD.getNodes()[0].getx();
		y1 = cD.getNodes()[0].gety();
		x2 = bD.getNodes()[i].getx(); // body node
		y2 = bD.getNodes()[i].gety();
		
		if(computeDistance2(x1,y1,x2,y2) < (overlapAllowanceD*radiusCD)*(overlapAllowanceD*radiusCD) )
		{
			fibroblasts[cellID].setMobility(cisD);
			//cout << "Judder!" << endl;
		}
	}	
}	

vector<double> TME::getMij()
{
	return Mij;
}

vector<FibreBox> TME::getFibreBoxes()
{
	return fibreBoxes;
}
std::vector<FibreBox> TME::getQuantBoxes()
{
	//find max Density
	//make map(int density, vector boxID)
	int max = 0; // to work out which is the maximum density of any fibreBox
	map<int, vector<FibreBox>> densityMap; // <density, boxesWithThatDensity>
	for(int i=0; i!= fibreBoxes.size(); i++)
	{
		if(woundSizeD!=0)
		{
			if(fibreBoxes[i].getWounded()==true)
			{
				int density = fibreBoxes[i].getDominantDensity();
				if(density > max)
				{
					max = density;
				}
				densityMap[density].push_back(fibreBoxes[i]);
			}
		}
		else
		{
			int density = fibreBoxes[i].getDominantDensity();
			if(density > max)
			{
				max = density;
			}
			densityMap[density].push_back(fibreBoxes[i]);
		}
	}

	// shuffle each vector<FibreBox>
	for(int i=0; i<=max; i++)
	{
		vector<FibreBox> f = densityMap[i];
		int entries = densityMap[i].size();
		if(entries>1) // only shuffle boxes for a given density if there is more than one box to shuffle
		{
			random_shuffle(f.begin(), f.end());
		}
		densityMap[i] = f;
	}
	
	//concatenate
	vector<FibreBox> fConcatenated;
	for(int i=max; i>=0; i--)
	{
		for(int j=0; j!=densityMap[i].size(); j++)
		{
			fConcatenated.push_back(densityMap[i][j]);
			if(fConcatenated.size() >= boxesToRecordD)
			{
				break;
			}
		}
		if(fConcatenated.size() >= boxesToRecordD)
		{
			break;
		}
	}
	
	return fConcatenated;
}

void TME::depositFibres(Fibroblast &cD)
{	
	// deposit matrix everywhere
	
	int previousBox = cD.getPreviousBox();
	fibreBoxes[previousBox].processNewFibres(cD.getTheta());
	
	// deposit matrix inside wound only
	/*double x1 = cD.getNodes()[0].getx();
	double y1 = cD.getNodes()[0].gety();
	
	if((x1-512)*(x1-512) + (y1-512)*(y1-512) < woundSizeD*woundSizeD) // if within 100000 of CCs then no good
	{
		int previousBox = cD.getPreviousBox();
		fibreBoxes[previousBox].processNewFibres(cD.getTheta());
	}*/	
}

double TME::sampleCIL()
{
	double changeInAngleDeg = 0; 
	double u1 = (double) rand() / (RAND_MAX); 
	int posVec = u1*CILdataHBD.size();
	changeInAngleDeg = CILdataHBD[posVec];	
	
	double changeInAngleRad = changeInAngleDeg*(M_PI/180);
	return changeInAngleRad;
}

vector<Fibroblast> TME::getFibroblasts()
{
	return fibroblasts;
}
void TME::readParameters(Parameter p)
{
	numberVoxBoxD = 1;
	rhoD = p.getRho();
	rhoCCD = 0;
	speedMeanD = p.getSpeedMean();
	speedSDD = p.getSpeedSD(); 
	polarityD = p.getPolarity();
	proliferationD = p.getProliferation(); 
	radiusCD = p.getRadiusC();
	radiusCCD = 1;
	CILaD = p.getCILa();
	w1D = p.getW1();
	w2D = p.getW2();
	w3D = 0;
	w4D = 0;
	w5D = 0;
	CILbD = p.getCILb();
	hetCILaD = 1;
	hetCILbD = 1;
	CILtypeD = p.getCILtype();
	lengthD = p.getLength();
	trajectoryFileNameD = p.getTrajectoryFileName();
	fibreGridLengthD =1;
	boxesToRecordD = 1;
	sourceXD = 1;
	sourceYD = 1;
	depRateD = 0;
	reRateD = 0;
	degRateD = 0;
	CILdataHBD = p.getCILdataHB();
	overlapAllowanceD = p.getOverlapAllowance();
	cisD = p.getCis();
	fakeRoseSDD = p.getFakeRoseSD();
	fakeRoseMeanD = p.getFakeRoseMean();
	IMTD = "N";
	IMID = 0;
	aspectRatioD = p.getAspectRatio();
	noBinsD = 1;
	headSizeD = p.getHeadSize();
	probCILD =p.getProbCIL();
	woundSizeD = 0;
}

void TME::readMatrix( int timeD)
{
	Mij.clear();
	string timeAsString = static_cast<ostringstream*>( &(ostringstream() << timeD) )->str();
	string fileName = "mCAT/mCAT_" + timeAsString + ".txt";
	ifstream matrixFile(fileName.c_str());
	string valueIJ;
	while(matrixFile >> valueIJ)
	{
		double mij = std::atof(valueIJ.c_str());
		Mij.push_back(mij);
	} 
}
