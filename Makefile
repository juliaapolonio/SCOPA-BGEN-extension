#SCOPA program

CC = g++

DEBUGFLAGS = -Wno-deprecated -O3 -lz

SCOPA:	main.cpp

	g++ ap.cpp chisquaredistr.cpp gammafunc.cpp global.cpp ialglib.cpp ibetaf.cpp igammaf.cpp main.cpp normaldistr.cpp regression.cpp structures.cpp studenttdistr.cpp tools.cpp $(DEBUGFLAGS) -o SCOPA
