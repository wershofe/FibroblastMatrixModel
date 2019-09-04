#ifndef NODE_H
#define NODE_H

class Node
{	
	public:
		Node();
		~Node();	
	
		void setx(double xD);
		double getx();
		void sety(double yD);
		double gety();
		void setNodeType(std::string nodeTypeD);
		std::string getNodeType();
				
	private:
	double x, y; 
	std::string nodeType; // Cells are made up of 4 nodes, each of which has a position and a type eg H, S, B, T		
};

#endif


