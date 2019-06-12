#C = cc -Wall
CC = g++
CFLAGS = -g -Wno-deprecated -fopenmp -std=c++11 
OBJS = Original1.o TME.o Fibroblast.o Contact.o Node.o Parameter.o FibreBox.o CC.o CAT.o
LIBS =  -lm -lstdc++ -lgomp
.SUFFIXES: .c .o
.c.o:
	$(CC) $(CFLAGS) -c $*.c -o $*.o
.cpp.o:
	$(CC) $(CFLAGS) -c $*.cpp -o $*.o
Original1: $(OBJS)
	$(CC) -o Original1 $(OBJS) $(LIBS)

clean:
	rm -f *.o
	
