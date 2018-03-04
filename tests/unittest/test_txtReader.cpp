#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "dEploidLasso.hpp"

class TestTextReader : public CppUnit::TestCase {
    CPPUNIT_TEST_SUITE( TestTextReader );
    CPPUNIT_TEST( testMainConstructor );
    CPPUNIT_TEST( testInitialization );
    CPPUNIT_TEST_SUITE_END();

  private:
    //double epsilon3;

  public:
    void setUp() {
        //this->epsilon3 = 0.000000000001;
    }

    void tearDown() {
    }

    void testMainConstructor(){

    }

    void testInitialization(){

    }
            //vector < vector <double> > matrix = TxtReader ("data/panel_chrom1.txt").matrix;
        //vector < double > wsaf = TxtReader ("data/PG0402-C_chrom1.wsaf").vec;
};

CPPUNIT_TEST_SUITE_REGISTRATION( TestTextReader );
