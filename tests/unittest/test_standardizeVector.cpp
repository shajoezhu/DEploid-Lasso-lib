#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "dEploidLasso.hpp"

class TestStandardizeVector : public CppUnit::TestCase {
    CPPUNIT_TEST_SUITE( TestStandardizeVector );
    CPPUNIT_TEST( testMainConstructor );
    CPPUNIT_TEST( testInitialization );
    CPPUNIT_TEST_SUITE_END();

  private:
        standardizeVector *dummy;

  public:
    void setUp() {
        dummy = new standardizeVector(TxtReader("data/myy.txt").vec);
    }

    void tearDown() {
        delete dummy;
    }

    void testMainConstructor(){
        //TxtReader panel("data/panel_chrom1.txt");
        //CPPUNIT_ASSERT_EQUAL(panel.matrix.size(), (size_t)625);
        //CPPUNIT_ASSERT_EQUAL(panel.matrix[0].size(), (size_t)94);
        //TxtReader wsaf("data/PG0402-C_chrom1.wsaf");
        //CPPUNIT_ASSERT_EQUAL(wsaf.vec.size(), (size_t)625);
    }

    void testInitialization(){
        //CPPUNIT_ASSERT_EQUAL(readerDummy_ptrX->matrix.size(), (size_t)10);
        //CPPUNIT_ASSERT_EQUAL(readerDummy_ptrX->matrix[0].size(), (size_t)2);
        //CPPUNIT_ASSERT_EQUAL(readerDummy_ptrY->vec.size(), (size_t)10);

    }
};

CPPUNIT_TEST_SUITE_REGISTRATION( TestStandardizeVector );
