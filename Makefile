CPP = g++
CFLAGS = -std=c++0x -Isrc/ -O3 -g -Wall -DNDEBUG
CFLAGS_dbg = -std=c++0x -Isrc/ -O3 -g -Wall


a_dbg.out: lasso_dbg.o lassoDBG_dbg.o main.cpp
	$(CPP) $(CFLAGS_dbg) -o a_dbg.out lasso_dbg.o lassoDBG_dbg.o main.cpp

lasso_dbg.o: src/dEploidLasso.cpp src/dEploidLasso.hpp
	$(CPP) $(CFLAGS_dbg) -c src/dEploidLasso.cpp -o lasso_dbg.o

lassoDBG_dbg.o: src/dEploidLasso.cpp src/dEploidLasso.hpp src/dEploidLassoDBG.cpp
	$(CPP) $(CFLAGS_dbg) -c src/dEploidLassoDBG.cpp -o lassoDBG_dbg.o


a.out: lasso.o lassoDBG.o main.cpp
	$(CPP) $(CFLAGS) -o a.out lasso.o lassoDBG.o main.cpp

lasso.o: src/dEploidLasso.cpp src/dEploidLasso.hpp
	$(CPP) $(CFLAGS) -c src/dEploidLasso.cpp -o lasso.o

lassoDBG.o: src/dEploidLasso.cpp src/dEploidLasso.hpp src/dEploidLassoDBG.cpp
	$(CPP) $(CFLAGS) -c src/dEploidLassoDBG.cpp -o lassoDBG.o




clean:
	rm a.out *.o
