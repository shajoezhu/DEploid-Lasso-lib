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

class TestDEploidLASSO : public CppUnit::TestCase {
    CPPUNIT_TEST_SUITE(TestDEploidLASSO);
    CPPUNIT_TEST(testValues);
    CPPUNIT_TEST(testValues2);
    CPPUNIT_TEST(testConstructor);
    CPPUNIT_TEST_SUITE_END();

 private:
        DEploidLASSO *dummy;
        double epsilon6;
        double epsilon5;
        double epsilon4;

 public:
    void setUp() {
        this->epsilon6 = 0.000001;
        this->epsilon5 = 0.00001;
        this->epsilon4 = 0.0001;
    }

    void tearDown() {
    }

    void testConstructor() {
        vector < vector <double> > matrix =
            lasso::TxtReader("data/myX.txt").matrix;
        vector < double > wsaf = lasso::TxtReader("data/myy.txt").vec;
        dummy = new DEploidLASSO(matrix, wsaf, 2);
        delete dummy;
        DEploidLASSO dummy0(matrix, wsaf, 0);
        DEploidLASSO dummy1(matrix, wsaf, 1);
        CPPUNIT_ASSERT_NO_THROW(dummy1.printResults());
    }



    void testValues() {
        vector < vector <double> > matrix =
            lasso::TxtReader("data/myX.txt").matrix;
        vector < double > wsaf = lasso::TxtReader("data/myy.txt").vec;

        DEploidLASSO dummy3(matrix, wsaf, 3);
    // Call:  glmnet(x = truex, y = y, lambda = 1/seq(3, 5, length.out = 3),
    // lower.limits = 0, type.gaussian = "naive")

    // Df   %Dev Lambda
    // [1,]  2 0.9949 0.3333
    // [2,]  2 0.9971 0.2500
    // [3,]  2 0.9982 0.2000
    // > myfit$beta
    // 2 x 3 sparse Matrix of class "dgCMatrix"
    //         s0       s1       s2
    // V1 9.446489 9.585377 9.668710
    // V2 4.439772 4.578661 4.661995
        vector <int> df({2, 2, 2});
        vector <double> dev({.9949, .9971, .9982});
        vector <double> lambda({.3333, .2500, .2000});
        vector <double> beta1({9.446489, 9.585377, 9.668710});
        vector <double> beta2({4.439772, 4.578661, 4.661995});
        for (size_t i = 0; i < 3; i++) {
            CPPUNIT_ASSERT_EQUAL(dummy3.df[i], df[i]);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(dummy3.devRatio[i], dev[i], epsilon4);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(dummy3.lambda[i], lambda[i], epsilon4);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(dummy3.beta[i][0], beta1[i], epsilon4);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(dummy3.beta[i][1], beta2[i], epsilon4);
        }
    }

    void testValues2(){
        vector < vector <double> > matrix =
            lasso::TxtReader("data/panel_chrom1.txt").matrix;
        vector < double > wsaf =
            lasso::TxtReader("data/PG0402-C_chrom1.wsaf").vec;

        DEploidLASSO dummy3(matrix, wsaf, 100);
    }
};

CPPUNIT_TEST_SUITE_REGISTRATION(TestDEploidLASSO);
