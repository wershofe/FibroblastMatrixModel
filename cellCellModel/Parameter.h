//28/10/16
#ifndef PARAMETER_H
#define PARAMETER_H
#include <vector>
#include <map>
#include <string>
#include <fstream>

class Parameter
{
	public:
		Parameter();
		~Parameter();

		void setAllParameters(std::string fileNameD);
		
		void readInCILdataHB(std::string CILdataHBnameD);
		void setNumberVoxBox();
		void setRho(int rhoD);
		void setSpeedMean(double speedMeanD);
		void setSpeedVar(double speedSDD);
		void setPolarity(double polarityD);
		void setProliferation(double proliferationD);
		void setRadiusC(double radiusCD);
		void setRadiusCC(double radiusCCD);
		void setCILa(double CILaD);
		void setHetCILa(double hetCILaD);
		void setW1(double w1D);
		void setW2(double w2D);
		void setW3(double w3D);
		void setW4(double w4D);
		void setW5(double w5D);
		void setCILb(double CILbD);
		void setHetCILb(double hetCILbD);
		void setCILtype(std::string CILtypeD);
		void setLength(int lengthD);
		void setTrajectoryFileName(std::string trajectoryFileNameD);
		void setFibreGridLength();
		void setSourceX(int sourceXD);
		void setSourceY(int sourceXD);
		int getNoBins();
		
		int getNumberVoxBox();
		int getRho();
		int getRhoCC();
		double getSpeedMean();
		double getSpeedSD();
		double getPolarity();
		double getProliferation();
		double getRadiusC();
		double getRadiusCC();
		double getCILa();
		double getHetCILa();
		double getW1();
		double getW2();
		double getW3();
		double getW4();
		double getW5();
		double getCILb();
		double getHetCILb();
		std::string getCILtype();
		int getLength();
		std::string getTrajectoryFileName();
		int getBoxesToRecord();
		int getFibreGridLength();
		int getSourceX();
		int getSourceY();
		int getDepRate();
		int getReRate();
		int getDegRate();
		std::string getIMT();
		int getIMI();
		std::vector<double> getCILdataHB();
		int getMaxFrame();
		double getOverlapAllowance();
		double getCis();
		double getFakeRoseSD();
		double getFakeRoseMean();
		int getAspectRatio();
		double getHeadSize();
		double getProbCIL();
		int getWoundSize();
		
	private:
		int numberVoxBox; // to do with chucking new cells on during proliferation 
		int rho; // initial number of cells
		int rhoCC; // initial number of CCs
		int maxFrame;
		double speedMean; // for now, each cell picks and fixes its speed from a Gaussian with given mean and speed
		double speedSD;
		double polarity; // intrinsic noise in cell motility
		double proliferation; // take the population doubling rate from data. 
		double radiusC; // radius of each node of the cell (apart from the head which is smaller);
		double radiusCC; // radius of CC - consists of only one node (the head)
		double CILa; // angle of collision at which CIL is activated
		double w1; // weighting of what cell i had been doing in previous timestep
		double w2; // weighting of neighbours (cell-cell adhesion);
		double w3; // weighting of ECM (currently 0);
		double w4; // heteroCG
		double w5; // ECM feedback
		double CILb; // the response angle when CIL is activated, 
		double hetCILa; // between Fibs and CCs
		double hetCILb;
		std::string CILtype; // the type of CIAb: fixed angle, drawn from uniform distribution~U[0,2*pi] ("N"); or from a Gaussian ~N(0,0.1);
		int length; // This is the size of the screen (normally 1024 pixels) and corresponds to the microscope zoom too.
		std::string trajectoryFileName;
		int boxesToRecord;
		int fibreGridLength;
		int sourceX;
		int sourceY;
		int depRate;
		int reRate;
		int degRate;
		std::string IMT; // initial matrix type
		int IMI; // initial matrix intensity
		std::vector<double> CILdataHB;
		double overlapAllowance;
		double cis; // change in speed upon VE trigger
		double fakeRoseSD;
		double fakeRoseMean;
		int aspectRatio;
		int noBins;
		double headSize;
		double probCIL;
		int woundSize;
};

#endif
