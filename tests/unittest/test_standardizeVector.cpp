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

class TestStandardizeVector : public CppUnit::TestCase {
    CPPUNIT_TEST_SUITE(TestStandardizeVector);
    CPPUNIT_TEST(testValues);
    CPPUNIT_TEST_SUITE_END();

 private:
        standardizeVector *dummy;
        double epsilon6;
        double epsilon5;

 public:
    void setUp() {
        dummy = new standardizeVector(lasso::TxtReader("data/myy.txt").vec);
        this->epsilon6 = 0.000001;
        this->epsilon5 = 0.00001;
    }

    void tearDown() {
        delete dummy;
    }

    void testValues(){
        // > normalized_y
        CPPUNIT_ASSERT_DOUBLES_EQUAL(dummy->mean, 7.50654, epsilon6);
        // var(y) * 9 /10
        CPPUNIT_ASSERT_DOUBLES_EQUAL(dummy->variance, 36.24513, epsilon5);
        // sqrt(var(y) * 9 /10)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(dummy->stdv, 6.020393, epsilon6);
        vector <double> normalized_y({0.3944360, -0.1315933, 0.3945837,
                                      -0.1312414, 0.3924621, -0.3943260,
                                      0.1313320, -0.3924229,
                                      0.1317995, -0.3950297});
        for (size_t i = 0; i < 10; i++) {
            CPPUNIT_ASSERT_DOUBLES_EQUAL(dummy->ret[i], normalized_y[i],
                                         epsilon6);
        }
    }
};

CPPUNIT_TEST_SUITE_REGISTRATION(TestStandardizeVector);
