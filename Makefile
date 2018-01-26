CPP = g++
CFLAGS = -std=c++0x -Isrc/ -O3 -g -Wall -DNDEBUG

a.out: lasso.o main.cpp
	$(CPP) $(CFLAGS) -o a.out lasso.o main.cpp

lasso.o: src/dEploidLasso.cpp src/dEploidLasso.hpp
	$(CPP) $(CFLAGS) -c src/dEploidLasso.cpp -o lasso.o
