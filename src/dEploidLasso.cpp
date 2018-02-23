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

#include "dEploidLasso.hpp"
#include <string>
#include <fstream>      // std::ifstream
#include <cmath>        // std::abs
#include <assert.h>     /* assert */
#include <limits> // std::numeric_limits<double>::infinity();


TxtReader::TxtReader(const char inchar[]){
    string fileName(inchar);
    ifstream in_file(inchar);
    string tmp_line;
    if ( in_file.good() ){
        getline ( in_file, tmp_line );
        while ( tmp_line.size() > 0 ){
            size_t field_start = 0;
            size_t field_end = 0;
            size_t field_index = 0;
            vector <double> contentRow;
            while ( field_end < tmp_line.size() ){
                field_end = min ( min ( tmp_line.find(',',field_start),
                                        tmp_line.find('\t',field_start) ),
                                  tmp_line.find('\n', field_start) );

                string tmp_str = tmp_line.substr( field_start, field_end - field_start );
                contentRow.push_back( strtod(tmp_str.c_str(), NULL) );
                if ( field_index == 0 ){
                    vec.push_back(strtod(tmp_str.c_str(), NULL));
                }
                field_start = field_end+1;
                field_index++;
            }
            this->matrix.push_back(contentRow);
            getline ( in_file, tmp_line );
        }
    }

    in_file.close();

}


standardizeVector::standardizeVector(vector <double> vec){
    //ym = mean(wsaf)
    //ytmp = (wsaf - ym)
    //ytmp / sqrt(sum(ytmp^2))
    //vector <double> ret;
    assert(ret.size() == 0);
    size_t nObs_ = vec.size();
    this->mean = sumOfVec(vec) / (double)nObs_;
    vector <double> mean_vec = vector <double> (nObs_, mean);
    vector <double> vec_diff = vecDiff(vec, mean_vec);
    vector <double> tmpProd = vecProd(vec_diff, vec_diff);
    this->stdv = sqrt(sumOfVec(tmpProd) / (double)nObs_);

    for ( double tmp : vec_diff ){
        ret.push_back(tmp/stdv);
    }
    assert(ret.size() == nObs_);
}


DEploidLASSO::DEploidLASSO(vector < vector <double> > &x, vector < double > &wsaf){
    this->nVars_ = x[0].size();
    this->nObs_ = x.size();

    // TODO
    // check for x size
    cout<< "Matrix size = "<< this->nObs_ << " " << this->nVars_ << endl;
    cout<< "Vector length = " << wsaf.size() << endl;

    // Initialize
    this->initialization();
    this->checkVariables(x);

    this->standarization(x, wsaf);

    this->productOfxy(); // TODO, needs work!!!!

    this->ix = vector <double> (nVars_, 0.0);

    this->nulldev_ = computeNullDev(x, wsaf);
    cout << "nulldev = " << this->nulldev_ << endl;
    for ( size_t i = 0; i < this->lambda.size(); i++){
        this->setLambdaCurrent(1.0 / (2.0+(double)i));
        this->lassoGivenLambda(this->ju, this->g, this->ix);

        // FETCH AND UPDATE THE CURRENT INFERENCE RESULTS
        this->beta.push_back(betaCurrent);
        this->lambda[i] = lambdaCurrent();
        this->intercept[i] = interceptCurrent();
        this->devRatio[i] = rsqCurrent();
        this->df[i] = dfCurrent();
    }
    cout<<"beta.size = "<< beta.size() <<endl;
}


DEploidLASSO::~DEploidLASSO(){

}


//DEploidLASSO::DEploidLASSO(){

//}

void DEploidLASSO::productOfxy(){

    //10780 do 10781 j=1,ni
            //if(ju(j).eq.0)goto 10781
            //g(j)=abs(dot_product(y,x(:,j)))
    //10781 continue

    this->g = vector <double> (nVars_, 0.0);

    for (size_t j = 0; j < this->nVars_; j++){
        if (this->ju[j] == 0){
            break;
        }
        //double tmpSum = 0;
        //for ( size_t i = 0; i < this->nObs_; i++ ){
            //tmpSum += standardized_y[i] * standardized_x[i][j];
        //}
        vector <double> tmp = vecProd(standardized_y, standardized_x_transposed[j]);
        this->g[j] = std::abs(sumOfVec(tmp));
    }

    assert(this->g.size() == nVars_);
}


void DEploidLASSO::checkVariables(vector < vector <double> > &x){
    //size_t nVars = x[0].size();
    //size_t nObs = x.size();

    this->ju = vector <double> (this->nVars_, 0.0);

//11060 do 11061 j=1,ni
            //ju(j)=0
            //t=x(1,j)
      //11070 do 11071 i=2,no
                  //if(x(i,j).eq.t)goto 11071
                  //ju(j)=1
                  //goto 11072
            //11071 continue
      //11072 continue
      //11061 continue
//11062 continue

    for (size_t i = 0; i < this->nVars_; i++){
        for (size_t ii = 1; ii < this->nObs_; ii++){
            if (x[ii][i] != x[0][i]){
                ju[i] = 1.0;
                break;
            }
        }
    }
}


void DEploidLASSO::standarization(vector < vector <double> > &x, vector < double > &y){
    // standarize x
    assert(standardized_x_transposed.size() == 0);
    assert(x_stdv.size() == 0);
    assert(x_mean.size() == 0);
    for (size_t i = 0; i < nVars_; i++){
        // Extract the ith variable
        vector <double> var_i;
        for ( size_t j = 0; j < this->nObs_; j++){
            var_i.push_back(x[j][i]);
        }
        standardizeVector vecX(var_i);
        x_stdv.push_back(vecX.stdv);
        x_mean.push_back(vecX.mean);
        standardized_x_transposed.push_back(vecX.ret);
    }
    assert(standardized_x_transposed.size() == nVars_);

    // standarize y
    standardizeVector vecY(y);
    standardized_y = vecY.ret;
    this->y_stdv = vecY.stdv;
    this->y_mean = vecY.mean;
}


void DEploidLASSO::initialization(size_t nLambda){
    this->lambda = vector <double> (nLambda, 0.0);
    this->intercept = vector <double> (nLambda, 0.0);
    this->devRatio = vector <double> (nLambda, 0.0);
    this->df = vector <int> (nLambda, 0);
    this->indexArray = vector <size_t> (nVars_, nVars_); // To initialize to nVars, as index are 0-based, index is from 0 to nVars-1

    this->lowerLimit = 0;
    this->upperLimit = std::numeric_limits<double>::infinity();

    this->maxIteration_ = 100000;
    this->thresh_ = 1e-7;
    this->absTol_ = 1e-4;
    this->dfmax_ = nVars_ + 1;
    this->npass_ = 0;
    this->mm = vector < size_t > (nVars_, (size_t)0);
    this->iz = 0;
    //this->nVars_ = nVars;
    //pmax = min(dfmax * 2+20, nvars)
    //this->nObs_ = nObs;

    //this->vp = vector <double> (nVars, 1.0);
    //this->vq = vector <double> (nVars, 1.0); // For now ...

  ////bta=beta
      ////omb=1.0-bta
      ////ix=0

      ////if(flmin .ge. 1.0)goto 10771
            ////eqs=max(eps,flmin)
            ////alf=eqs**(1.0/(nlam-1))
////10771 continue

      ////rsq=0.0
      ////a=0.0
      ////mm=0
      ////nlp=0
      ////nin=nlp
      ////iz=0
      ////mnl=min(mnlam,nlam)
      ////alm=0.0


}


double DEploidLASSO::computeNullDev(vector < vector <double> > &x, vector < double > &wsaf){
    double ybar = sumOfVec(wsaf) / (double)wsaf.size();
    vector <double> ybar_vec = vector <double> (wsaf.size(), ybar);
    vector <double> diff = vecDiff(wsaf, ybar_vec);
    vector <double> tmpSq = vecProd(diff, diff);
    return (double) sumOfVec( tmpSq );
}


void DEploidLASSO::lassoGivenLambda(vector < double > ju, vector <double> g, vector <double> &ix){
    /* USE THE FOLLOWING VARIABLES
     *
     * lambdaCurrent()
     * standardized_x_transposed
     * standardized_y
     *
     */
    //cout << lambda << endl;
    assert(nVars_ > 0);
    assert(nObs_ > 0);

    /*
     * INITIALIZATION
     */
    this->betaCurrent = vector <double> (this->nVars_, 0.0);
    this->coefficentCurrent = vector <double> (this->nVars_, 0.0);
    this->setDfCurrent(0);
    this->setRsqCurrent(0.0);
    this->setInterceptCurrent(0.0);

    //vector <double> vp(this->nVars_, 1.0);
    this->jz = 1;

    // ulam is user defined lambdas...
    // no = number of observations
    // ni = number of predictor variables

    //c   parm = penalty member index (0 <= parm <= 1)
    //c        = 0.0 => ridge
    //c        = 1.0 => lasso
    double beta = 1.0;
    double bta=beta; // lasso part
    double omb=1.0-bta; // ridge part, omb = 0
    double rsq=0.0; //TODO, just replace as rsq0

    // Initializations
    double alm0 = 0.0;
    //double alm  = lambda;
    double flmin = 1.0; // this is defined in glmnet.R

    //dem=alm*omb
    //ab=alm*bta
    //rsq0=rsq
    //jz=1
    //tlam=bta*(2.0*alm-alm0)

    //double dem = lambdaCurrent() * omb;
    double dem = 0;
    double ab = lambdaCurrent() * bta;
    double rsq0 = rsq;
    this->jz=1;
    double tlam=bta*(2.0 * lambdaCurrent() - alm0);


      //rsq=0.0
      //a=0.0
      //mm=0
      //nlp=0
      //nin=nlp
      //iz=0
      //mnl=min(mnlam,nlam)
      //alm=0.0

    //10860 do 10861 k=1,ni
                //if(ix(k).eq.1)goto 10861
                //if(ju(k).eq.0)goto 10861
                //if(g(k).gt.tlam*vp(k)) ix(k)=1
          //10861 continue
    //10862 continue

    for ( size_t k = 0; k < this->nVars_; k++){
        if (ix[k] == 1) {break;}
        if (ju[k] == 0) {break;}
        if (g[k] > (tlam*1.0)){
            ix[k] = 1.0;
        }
    }


this->updateWithNewVariables(ab);
                  //10870 continue
                  //10871 continue
                        //if(iz*jz.ne.0) go to 10360
                  //10880 continue
          //10992 continue
                //jz=0
                //goto 10871
          //10872 continue


    // Rescale coefficients
    this->rescaleCoefficents();

    // Compute intercept
    this->computeIntercept();

    // Map coefficients, coefficient -> beta
    this->coefficentToBeta();
}


void DEploidLASSO::updateWithNewVariables(double ab){
    double ixx;
    while ( true ){
        if ( iz*jz != 0 ){
            iz = 1;
            updateWithTheSameVariables(ab);
        } else {
                  this->npass_++;
                  double dlx=0.0;
                  for (size_t k = 0; k < nVars_; k++){
                      if (ix[k] == 0){ break; }
                      double ak = this->coefficentCurrent[k];
                        double gk = this->updateAkReturnGk(k, ak, ab);
                      if(coefficentCurrent[k] == ak) {
                          break;
                      }

                      if (mm[k] == 0){
                          nin += 1;
                          if (nin > nVars_){
                              break;
                          }
                          mm[k] = nin;
                          indexArray[nin] = k;
                      }
                      this->updateY(k, gk, ak, dlx);
                  }

for ( size_t k = 0; k < nVars_; k++){
    if (ix[k] == 1){break;}
    if (ju[k] == 0){break;}

    vector <double> tmp = vecProd(standardized_y, standardized_x_transposed[k]);
g[k] = std::abs(sumOfVec(tmp));
    if (g[k] <= ab){break;}
    ix[k] = 1;
    ixx = 1;

    }

                   //10940 do 10941 k=1,ni
                        //if(ix(k).eq.1)goto 10941
                        //if(ju(k).eq.0)goto 10941
                        //g(k)=abs(dot_product(y,x(:,k)))
                        //if(g(k) .le. ab*vp(k))goto 10961
                        //ix(k)=1
                        //ixx=1
                  //10961 continue
                  //10941 continue
                  //10942 continue

                        if(ixx != 1) {break;}
            if (npass_ > maxIteration_){
                return;
            }
        }
    }

}



void DEploidLASSO::updateWithTheSameVariables(double ab){
    while ( npass_ < this->maxIteration_ ){
        this->npass_++;
        double dlx = 0.0;
        for ( size_t l = 0; l < this->nin; l++){
            size_t k = indexArray[l];
            double ak = this->coefficentCurrent[k];
            double gk = this->updateAkReturnGk(k, ak, ab);
            if(coefficentCurrent[k] == ak) {
            break;
            }
            dlx = this->updateY(k, gk, ak, dlx);
        }
        if (dlx < thresh_){
            break;
        }
    }
    jz = 0;
}


double DEploidLASSO::updateAkReturnGk(size_t k, double ak, double ab){
    cout<<"here"<<endl;
    vector <double> gk_vec = vecProd(this->standardized_y,
                                     this->standardized_x_transposed[k]);
    double gk = sumOfVec(gk_vec);
    double u = gk+ak;
    double u_sign = 1.0;

    if ( u < 0 ){
        u_sign = -1.0;
    }
    double v=abs(u)-1*ab;
    coefficentCurrent[k] = 0.0;
    if(v > 0.0){
      coefficentCurrent[k] = max(lowerLimit,
                 min(upperLimit,
                    v*u_sign));
    }
    return gk;
}


double DEploidLASSO::updateY(size_t k, double gk, double ak, double dlx){
    double del = this->coefficentCurrent[k] - ak;
    this->rsqCurrent_ += del * (2.0*gk-del);
    for ( size_t i = 0; i < this->nObs_; i++ ){
        standardized_y[i] -= del * standardized_x_transposed[k][i];
    }

    return max(del*del, dlx);
}


void DEploidLASSO::computeIntercept(){
    this->setInterceptCurrent(0.0);

    if (ninCurrent_ == 0){
        return;
    }

    //a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))

    double y_remaining = this->y_mean;
    for ( size_t i = 0; i < (size_t)this->ninCurrent_; i++ ){
        y_remaining -= coefficentCurrent[i] * x_mean[indexArray[i]];
    }
    this->setInterceptCurrent(y_remaining);
}


void DEploidLASSO::rescaleCoefficents(){
    for ( size_t i = 0; i < (size_t)this->ninCurrent_; i++ ){
        this->coefficentCurrent[i] *= y_stdv;
        this->coefficentCurrent[i] /= x_stdv[indexArray[i]];
    }

    ////10620 do 10621 k=1,lmu
            ////alm(k)=ys*alm(k)
            ////nk=nin(k)
      //10630 do 10631 l=1,nk
                  //ca(l,k)=ys*ca(l,k)/xs(ia(l))
            //10631 continue
      //10632 continue
}


void DEploidLASSO::coefficentToBeta(){

}
