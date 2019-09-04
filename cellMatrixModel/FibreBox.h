//13/6/16
#ifndef FIBREBOX_H
#define FIBREBOX_H
#include <vector>
class FibreBox
{	
	public:
		FibreBox(int xD, int yD, int boxD, int dominantIndexD, int dominantDensityD, int depRateD, int reRateD, int degRateD, int noBinsD, bool woundedD);
		FibreBox(int xD, int yD, int boxD, std::vector<int> densities, int depRateD, int reRateD, int degRateD, int noBinsD, bool woundedD);
		~FibreBox();	
		
		void setX(int xD);
		void setY(int yD);
		void setBox(int boxD);
		
		int getX();
		int getY();
		int getBox();
		
		void setBins(std::vector<double> binsD);
		void setDensityBins(std::vector<int> densityBinsD);
		std::vector<double> getBins();
		std::vector<int> getDensityBins();
		
		void processNewFibres(double thetaD);
			int computeBin(double thetaD);
			void degradeAll();
			void rearrangeNeighbours(int indexD);
			void depositFibre(int indexD);
			void setDominantIndex();
		
		
		void degradeBin(int i);
		void rearrangeBin(int i);
		
		
		void setDominantDensity(int dominantDensityD);
		void setDominantOrientation(double dominantOrientationD);
		
		int getDominantIndex();
		int getDominantDensity();
		double getDominantOrientation();
		bool getWounded();
						
	private:
	int x, y, box; // which box is it
	int dominantIndex; // index of dominant bin (to get orientation and density)
	int dominantDensity;
	double dominantOrientation;
	std::vector<double> bins;
	std::vector<int> densityBins; 	
	bool wounded;	
	
	int depositionRate;
	int rearrangementRate;
	int degradationRate;
	int noBins;
};

#endif

