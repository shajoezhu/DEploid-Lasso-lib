/*
 * dEploid is used for deconvoluting Plasmodium falciparum genome from
 * mix-infected patient sample. DEploid-Lasso-lib is a submodule for
 * choosing the appropriate reference panel using the LASSO algorithm.
 *
 * Copyright (C) 2018 University of Oxford
 *
 * Author: Sha (Joe) Zhu
 *
 * This file is part of DEploid-Lasso-lib.
 *
 * dEploid is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "dEploidLasso.hpp"

class TestTextReader : public CppUnit::TestCase {
    CPPUNIT_TEST_SUITE( TestTextReader );
    CPPUNIT_TEST( testMainConstructor );
    CPPUNIT_TEST( testInitialization );
    CPPUNIT_TEST_SUITE_END();

  private:
        TxtReader *readerDummy_ptrX;
        TxtReader *readerDummy_ptrY;

  public:
    void setUp() {
        readerDummy_ptrX = new TxtReader("data/myX.txt");
        readerDummy_ptrY = new TxtReader("data/myy.txt");
    }

    void tearDown() {
        delete readerDummy_ptrX, readerDummy_ptrY;
    }

    void testMainConstructor(){
        TxtReader panel("data/panel_chrom1.txt");
        CPPUNIT_ASSERT_EQUAL(panel.matrix.size(), (size_t)625);
        CPPUNIT_ASSERT_EQUAL(panel.matrix[0].size(), (size_t)94);
        TxtReader wsaf("data/PG0402-C_chrom1.wsaf");
        CPPUNIT_ASSERT_EQUAL(wsaf.vec.size(), (size_t)625);
    }

    void testInitialization(){
        CPPUNIT_ASSERT_EQUAL(readerDummy_ptrX->matrix.size(), (size_t)10);
        CPPUNIT_ASSERT_EQUAL(readerDummy_ptrX->matrix[0].size(), (size_t)2);
        CPPUNIT_ASSERT_EQUAL(readerDummy_ptrY->vec.size(), (size_t)10);
    }
};

CPPUNIT_TEST_SUITE_REGISTRATION( TestTextReader );
