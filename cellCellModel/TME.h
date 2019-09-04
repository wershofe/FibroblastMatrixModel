// 7/11/17
//---------------------------------
/*
w1 = cell persistence
w2 = contact guidance between CAFs
w3 = CAT
w4 = interaction with CCs
w5 = CAF-ECM feedback
*/
//---------------------------------
#ifndef TME_H
#define TME_H
#include <vector>
#include <map>
#include "Fibroblast.h"
#include "Contact.h"
#include "Node.h"
#include <string>
#include "Parameter.h"
#include "FibreBox.h"
#include "CC.h"
#include <omp.h>
class TME
{
	public:
		TME();
		~TME();
		
		double computeDistance2(Fibroblast cD, Fibroblast bD); // computes distance including periodic boundary conditions
		double computeDistance2(double x1, double y1, double x2, double y2); // between specified pair of coordinates
		void boundaryChecks(double &x1D, double &x2D, double &x3D, double &x4D); // apply periodic boundaries in modulo sense for distances
		
		
		void initialiseTME();
			void initialiseFibroblasts();
			void initialiseFibreBoxes();
				bool woundMatrix(int &x, int &y);
			void initialiseCCs();
		
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		void updateTME(int &timeD); // NB different to update in Original1.cpp, which concerns openGL
		
			std::vector<double> assembleNewAngles(int &time); // makes vector of new angles - simultaneous
				//step 1
				std::vector<Contact> findContacts(Fibroblast &cD, int &timeD);
					void checkContacts(Fibroblast &cD, Fibroblast &bD, int iD, int jD, int &timeD, std::vector<Contact> &contactsD); //quick voxel
						void checkCIL(Fibroblast &cD, Fibroblast &bD, int iD, int jD, Contact conD, int &timeD); // check if contact triggers CIL
							double drawFromRose(); //CIL response either fixed or drawn from a prob distribution
				//step 2
				std::vector<Contact> findContactsCC(Fibroblast &cD, int &timeD);
					void checkContactsCC(Fibroblast &cD, CC  &ccD, int iD, int &timeD, std::vector<Contact> &contactsCCD);
						void checkHetCIL(Fibroblast &cD, CC &bD, int iD, Contact conD, int &timeD);
				
				///////////////////////////////////////////////////////////////////////////////////////////////////////		
				//step 3
				double computeAngle(Fibroblast &cD, std::vector<Contact> &contactsD,  std::vector<Contact> &contactsCCD, int &timeD); 

					// All these functions go inside computeAngle function
					//step a
					void computeW1Component(Fibroblast &cD, std::vector<Contact> &contactsD, double &w1ComponentD, int &timeD);
					void computeW2Component(Fibroblast &cD, std::vector<Contact> &contactsD, int &timeD);
						void modifyW2Element(Contact &contactD, double &w2ElementD);//13-2-17 so w2 is less than pi/2 from theta_i 
						void checkHeadHeadcontacts(Contact &contactD, double &w2ElementD, Fibroblast &cD);//nematic
						void checkOtherContacts(Contact &contactD, double &w2ElementD, Fibroblast &cD);
						void reverseAngle(Contact &contactD, double &w2ElementD); // nematic (also used in computeW4Component)
					void computeW3Component(Fibroblast &cD, int &timeD); //CAT
					void computeW4Component(Fibroblast &cD, std::vector<Contact> &contactsCCD, int &timeD); //CCs
						void modifyW4Element(Contact &contactD, double &w4ElementD);
						void checkHetNematic(Contact &contactCCD, double &w4ElementD, double &lowerBoundD, double &upperBoundD);
					void computeW5Component(Fibroblast &cD, int &timeD, int &noConsideredFibres, double &w5ComponentD);//CAF-ECM feedback
					//step b
					double determineNewAngle(Fibroblast &cD, double &w1ComponentD, double &w5ComponentD, int &noConsideredFibresD); //main model
					
					
			void blockFromCCs();
			void judderOverlaps();
			/*
			void stopOverlaps();
				void findNucAllContacts(); // both nucHead and nucNuc contacts
					void checkNucNucOverlap(Fibroblast &cD, Fibroblast &bD); //inside findNucAllContacts
					void checkNucHeadOverlap(Fibroblast &cD, Fibroblast &bD); //inside findNucAllContacts
			*/
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				
		std::vector<Fibroblast> getFibroblasts();
		std::vector<int> getFreeBox(); // proliferation - when need a new fibroblast, pick fibroblast at random & chuck new fibroblast in adjacent freebox
		int getNumberVoxBox(); //psuedo-grid into which new daughter fibroblasts can appear
		void readParameters(Parameter p);
		
		std::vector<FibreBox> getFibreBoxes(); 
		std::vector<FibreBox> getQuantBoxes(); 
		void depositFibres(Fibroblast &cD); 
				
		std::vector<CC> getCCs();
		
		void readMatrix( int timeD);
		double computeLevelCAT(Fibroblast &cD);
		std::vector<double> getMij();
		
		void setBoxToFibs(); // assigns each fibroblasts to a voxBox
		void setBoxToVoxFibs(); // rewrites map between voxBox number and all the fibroblasts in that box and adjacent boxes
			void computeAdjacentBoxes(int i); // want to know all fibroblasts in a voxel box and adjacent boxes - only those fibroblasts can be touching
		
		double sampleCIL();
		void addFibroblast(int &timeD);
		
		void freezeCell(int numCell);
		void checkVE(Fibroblast &cD, Fibroblast &bD);
	private: 
		std::vector<Fibroblast> fibroblasts, proposals; //proposals gives candidate new fibroblast positions - then check for overlaps and fix these for "fibroblasts".
		int numberVoxBox;
		double orderParameter, avgPersistence; // not using either right now
		
		int gridBoxLengthD; // number of gridboxes wanted for ECM deposition (total number is gridBoxLength squared); - not using this at the moment
		int numberVoxBoxD; // to do with chucking new fibroblasts on during proliferation 
		int rhoD; // initial number of fibroblasts
		int rhoCCD; // initial number of CCs
		double speedMeanD; // for now, each fibroblast picks and fixes its speed from a Gaussian with given mean and speed
		double speedSDD;
		double polarityD; // intrinsic polarity in fibroblast motility
		double proliferationD; // take the population doubling rate from data. 
		double radiusCD; // radius of each node of the fibroblast (apart from the head which is smaller);
		double radiusCCD; // radius of CC, which consists of just the head node
		double CILaD; // angle of collision at which CIL is activated
		double w1D; // weighting of what fibroblast i had been doing in previous timestep
		double w2D; // weighting of neighbours (fibroblast-fibroblast adhesion);
		double w3D; // CAT
		double w4D; // interaction between CCs and fibs
		double w5D; // CAF-ECM feedback
		double CILbD; // the response angle when CIL is activated, 
		double hetCILaD; // heterotypic angle between CAFs and CCs at which CIL is triggered
		double hetCILbD;
		
		std::string CILtypeD; // the type of CIAb: fixed angle, drawn from uniform distribution~U[0,2*pi] ("N"); or from a Gaussian ~N(0,0.1);
		int lengthD; // This is the size of the screen (normally 1024 pixels) and corresponds to the microscope zoom too.
		std::string trajectoryFileNameD;
		int frame;
		double jamRatio;
		std::vector<double> CILdataHBD; //sampling from roseplots
		double overlapAllowanceD;
		double cisD;
		double fakeRoseSDD;
		double fakeRoseMeanD;
		int aspectRatioD; // number of body nodes
		double headSizeD;
		double probCILD;
		int woundSizeD;
		
		std::vector<Contact> nucHeadContacts; // second level of volume exclusion (first being CIL)
		std::vector<Contact> nucNucContacts; // third level of volume exclusion
		
		std::vector<FibreBox> fibreBoxes; // ECM is deposited in boxes ie discretised to grid
		std::vector<FibreBox> quantBoxes; // ECM is deposited in boxes ie discretised to grid
		int boxesToRecordD;
		int fibreGridLengthD; // number of boxes (which determines box size)
		int sourceXD; // source of CAT, which is also for now the centre of the CCs
		int sourceYD;
		std::vector<CC> CCs;
		
		std::vector<double> Mij; // this is a vector storing the matrix of CAT levels
		
		std::map<int, std::vector<Fibroblast> > boxToFibs; // for each box, lists fibroblasts in that box - for voxel
		std::map<int, std::vector<Fibroblast> > boxToVoxFibs; // for each box, lists fibroblasts in that box and in all adjacent boxes
		
		int depRateD; //ECM deposition, rearrangement and degradation
		int reRateD;
		int degRateD;
		std::string IMTD; // initial matrix type (aligning or non-aligning)
		int IMID; // initial matrix intensity
		int noBinsD;
};

#endif
	
