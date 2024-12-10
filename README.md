# FibroblastMatrixModel

Folder cellMatrixModel pertains to the paper Wershof et al. 2019.

Folder cellCellModel pertains to the paper Park et al. 2019.

Contact.cpp = Controls nature of collisions between fibroblasts
CAT.cpp = Optional chemoattractant
CC.cpp = Optional cancer cells.
Fibrebox.cpp = matrix grid
Node.cpp = Fibroblasts are made up of nodes
TME.cpp = Tumour microenvironment
Original1.cpp = Controlling output files. Launch ./Original1 to run code
parametersToRead.txt = altering input parameters


Compile with
g++ -o Original1 CAT.cpp CC.cpp Contact.cpp FibreBox.cpp Fibroblast.cpp Node.cpp Original1.cpp Parameter.cpp TME.cpp

TME.cpp:
void TME::initialiseFibreBoxes()
FibreBox(int xD, int yD, int boxD, int dominantIndexD, int dominantDensityD, int depRateD, int reRateD, int degRateD, int noBinsD, bool woundedD)
IMT == ‘A’ | ’N’ | filename
If IMT == ‘A’ | ’N’ then IMI is referenced (to figure out how much initial synthetic matrix to put down)


FibreBox.cpp
FibreBox(int xD, int yD, int boxD, int dominantIndexD, int dominantDensityD, int depRateD, int reRateD, int degRateD, int noBinsD, bool woundedD)

[FILENAME]_fibrebox.txt Records:
box, frame, x,y, orientation, density
