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

template <typename T>
vector <T> matrixTimesVec( vector < vector < T > > &x, vector <double> &b ){
    vector <double> ret (x.size(), 0.0);
    for ( size_t i = 0; i < x.size(); i++ ){
        for ( size_t k = 0; k < x[i].size(); k++){
            ret[i] += x[i][k] * b[k];
        }
    }
    return ret;
}


template <typename T>
vector <T> vecDiff ( vector<T> &vecA, vector<T> &vecB ){
    //assert(vecA.size() == vecB.size());
    vector <T> difference (vecA.size(), (T)0);
    for ( size_t i = 0; i < vecA.size(); i++ ){
        difference[i] = vecA[i] - vecB[i];
    }
    return difference;
}


template <typename T>
vector <T> vecSum ( vector<T> &vecA, vector<T> &vecB ){
    assert(vecA.size() == vecB.size());
    vector <T> tmpSum (vecA.size(), (T)0);
    for ( size_t i = 0; i < vecA.size(); i++ ){
        tmpSum[i] = vecA[i] + vecB[i];
    }
    return tmpSum;
}


template <typename T>
vector <T> vecProd ( vector<T> &vecA, vector<T> &vecB ){
    //assert(vecA.size() == vecB.size());
    vector <T> tmpProd (vecA.size(), (T)0);
    for ( size_t i = 0; i < vecA.size(); i++ ){
        tmpProd[i] = vecA[i] * vecB[i];
    }
    return tmpProd;
}


template <typename T>
T sumOfVec( vector <T>& array ){
    T tmp = 0;
    for (auto const& value: array){
        tmp += value;
    }
    return tmp;
}


double computeNullDev(vector < vector <double> > &x, vector < double > &wsaf);


class LASSOgivenLambda{
  friend class DEploidLASSO;
  private:
    int maxIteration_;
    double thresh_;
    double absTol_;

    size_t nObs_;
    size_t nVars_;
    void initialization(size_t nObs);
    vector <double> beta;
    //double lambda; // size of nLambda
    double devRatio;
    double intercept;
    int df;
    vector <double> vp; // penalty
    vector <double> vq;
    LASSOgivenLambda(vector < vector <double> > &x, vector < double > &wsaf, double lambda);
    ~LASSOgivenLambda(){}
};



class DEploidLASSO{
#ifdef UNITTEST
  friend class TestDEploidLASSO;
#endif
  public:
    //DEploidLASSO();
    DEploidLASSO(vector < vector <double> > &x, // nObs x nVariable
                 vector < double > &wsaf);
    ~DEploidLASSO();

  private:
    void initialization(size_t nLambda = 100);
    vector < vector <double> > beta;
    vector < double > lambda; // size of nLambda
    vector < double > devRatio;
    vector < double > intercept;
    vector < int > df;
    double nulldev_;

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
