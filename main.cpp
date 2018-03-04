#include <iostream> // std::cout
#include "dEploidLasso.hpp"

using namespace std;

int main(int argc, char *argv[]){
    try {
        //vector < vector <double> > matrix = TxtReader ("data/panel_chrom1.txt").matrix;
        //vector < double > wsaf = TxtReader ("data/PG0402-C_chrom1.wsaf").vec;
        vector < vector <double> > matrix = TxtReader ("data/myX.txt").matrix;
        vector < double > wsaf = TxtReader ("data/myy.txt").vec;
        DEploidLASSO dummy(matrix, wsaf, 1);
        dummy.printResults();
        matrix.clear();
        wsaf.clear();
        return EXIT_SUCCESS;
    }
    catch (const exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
}
