#include <iostream> // std::cout
#include "dEploidLasso.hpp"

using namespace std;

int main(int argc, char *argv[]){
    try {
        //vector < vector <double> > matrix = TxtReader ("tests/panel_chrom1.txt").matrix;
        //vector < double > wsaf = TxtReader ("tests/PG0402-C_chrom1.wsaf").vec;
        vector < vector <double> > matrix = TxtReader ("tests/myX.txt").matrix;
        vector < double > wsaf = TxtReader ("tests/myy.txt").vec;
        DEploidLASSO dummy(matrix, wsaf);
        return EXIT_SUCCESS;
    }
    catch (const exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
}
