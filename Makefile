#SCOPA program

CC = g++

DEBUGFLAGS = -Wno-deprecated -O3 -lz

ALGLIB = $(wildcard ALGLIB/*.cpp)
TCLAP = $(wildcard TCLAP/*.cpp)
TOOLS = $(wildcard TOOLS/*.cpp)

SCOPA:	main.cpp

	g++ $(ALGLIB) $(TCLAP) global.cpp main.cpp $(TOOLS) $(DEBUGFLAGS) -o SCOPA
