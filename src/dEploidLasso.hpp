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
#include <sstream>
//#include <exceptions.hpp>

#ifndef NDEBUG
#define dout std::cout << "   "
#else
#pragma GCC diagnostic ignored "-Wunused-value"
#define dout 0 && std::cout
#endif

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


//class LASSOgivenLambda{
  //friend class DEploidLASSO;
  //private:
    //int maxIteration_;
    //double thresh_;
    //double absTol_;

    //size_t nObs_;
    //size_t nVars_;
    //size_t dfmax_;
    //void initialization(size_t nObs, size_t nVars);
    //vector <double> beta;
    ////double lambda; // size of nLambda
    //double devRatio;
    //double intercept;
    //int df;
    //vector <double> vp; // penalty
    //vector <double> vq;
    //LASSOgivenLambda(vector < vector <double> > &x, vector < double > &y, double lambda, vector < double > ju, vector <double> g, vector <double> &ix);
    //~LASSOgivenLambda(){}
//};


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


struct standardizeVector {
    vector <double> ret;
    double mean;
    double stdv;
    double variance;

    standardizeVector(vector <double> vec);
};


class DEploidLASSO{
#ifdef UNITTEST
  friend class TestDEploidLASSO;
#endif
  public:
    //DEploidLASSO();
    DEploidLASSO(vector < vector <double> > &x, // nObs x nVariable
                 vector < double > &y);
    ~DEploidLASSO();

  private:

  // FUNCTIONS
   // COMMON
    void initialization(size_t nLambda = 3);
    void standarization(vector < vector <double> > &x, vector < double > &y);
    void checkVariables(vector < vector <double> > &x);
    void productOfxy();
    double computeNullDev(vector < vector <double> > &x, vector < double > &y);

   // FOR EACH LAMBDA UPDATE
    void lassoGivenLambda();

  // VARIABLES, GETTERS AND SETTERS
   // COMMON
    size_t nObs_;
    size_t nVars_;
    vector < vector <double> > standardized_x_transposed; // nVariable x nObs
    vector <double> standardized_y;
    vector < vector <double> > beta; // nLambda x nVars
    vector < double > lambda; // size of nLambda
    vector < size_t > indexArray;
    vector < size_t > mm; // indicator, that kth variable is already in use
    size_t nin; // number of variables in use
    vector < double > devRatio;
    vector < double > intercept;
    int maxIteration_;
    double thresh_;
    //double absTol_;
    vector < int > df;
    double nulldev_;
    vector <double> ju;
    vector <double> g;
    vector <double> x_mean;
    vector <double> x_stdv;
    vector <double> x_variance;
    double y_stdv;
    double y_mean;
    size_t dfmax_;
    vector <double> ix;
    int npass_;

    double lowerLimit;
    double upperLimit;

   // FOR EACH LAMBDA UPDATE
   // VARIABLES
    vector < double > betaCurrent; // size of nVars
    vector < double > coefficentCurrent; // size of nVars

    double lambdaCurrent_;
    void setLambdaCurrent ( const double setTo ){ this->lambdaCurrent_ = setTo; }
    double lambdaCurrent() const { return this->lambdaCurrent_; }

    double lambdaPrevious_;
    void setLambdaPrevious ( const double setTo ){ this->lambdaPrevious_ = setTo; }
    double lambdaPrevious() const { return this->lambdaPrevious_; }

    double lambdaCurrentScaled_;
    void setLambdaCurrentScaled ( const double setTo ){ this->lambdaCurrentScaled_ = setTo; }
    double lambdaCurrentScaled() const { return this->lambdaCurrentScaled_; }

    double rsqCurrent_;
    void setRsqCurrent ( const double setTo ){ this->rsqCurrent_ = setTo; }
    double rsqCurrent() const { return this->rsqCurrent_; }

    double interceptCurrent_;
    void setInterceptCurrent ( const double setTo ){ this->interceptCurrent_ = setTo; }
    double interceptCurrent() const { return this->interceptCurrent_; }

    int dfCurrent_;
    void setDfCurrent ( const double setTo ){ this->dfCurrent_ = setTo; }
    double dfCurrent() const { return this->dfCurrent_; }

    int ninCurrent_;
    int iz, jz;
   // FUNCTIONS
    void computeIntercept();
    void rescaleCoefficents();
    void coefficentToBeta();
    void updateCoefficient(size_t k, double previousCoefficentValue, double gk);
    double updateYReturnDel(size_t k, double gk, double ak);
    void updateWithNewVariables();
    void updateWithTheSameVariables();
    void updatingCore();
    double computeGk(vector<double> &y, vector<double> &x);
    double computeGk_abs(vector<double> &y, vector<double> &x);
    void chooseVariables(double tlam);

   // Debug tools
    bool print_normalized_struff();
    bool print_initial_gk();
    bool print_homogeneous_input();
};


#endif
