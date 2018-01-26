/*
 * dEploid is used for deconvoluting Plasmodium falciparum genome from
 * mix-infected patient sample.
 *
 * Copyright (C) 2016-2017 University of Oxford
 *
 * Author: Sha (Joe) Zhu
 *
 * This file is part of dEploid.
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

#include <vector>
#include <iostream>
//#include <exceptions.hpp>
#include <sstream>

#ifndef LASSO
#define LASSO

using namespace std;

class DEploidLASSO{
#ifdef UNITTEST
  friend class TestDEploidLASSO;
#endif
  public:
    DEploidLASSO();
    DEploidLASSO(vector < vector <double> > &x, // nObs x nVariable
                 vector < double > &wsaf);
    ~DEploidLASSO();

  private:
    int maxIteration_;
    double realTol_;
    double absTol_;


    void initialization();
    size_t nObs_;
    vector < double > lambda; // size of
    vector < vector <double> > beta;
    vector <double> currentBeta;
    vector < double > devRatio;
    vector < int > df;
};


class TxtReader{
  friend class DEploidLASSO;

  public:
    vector < vector < double > > matrix;
    vector < double > vec;
    TxtReader (const char inchar[]);
    //virtual void readFromFile( const char inchar[] ){ this->readFromFileBase( inchar ); };
    //void readFromFileBase( const char inchar[] );
    ~TxtReader(){ };

  private:
};

#endif
