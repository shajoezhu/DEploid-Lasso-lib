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

double computeNullDev(vector < vector <double> > &x, vector < double > &wsaf){
    double ybar = sumOfVec(wsaf) / (double)wsaf.size();
    vector <double> ybar_vec = vector <double> (wsaf.size(), ybar);
    vector <double> diff = vecDiff(wsaf, ybar_vec);
    vector <double> tmpSq = vecProd(diff, diff);
    return (double)sumOfVec( tmpSq );
}


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


DEploidLASSO::~DEploidLASSO(){

}


//DEploidLASSO::DEploidLASSO(){

//}


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
    // need to standarise wsaf first
    this->productOfxy(x, wsaf); // TODO, needs work!!!!

    this->ix = vector <double> (nVars_, 0.0);

    this->nulldev_ = computeNullDev(x, wsaf);
    cout << "nulldev = " << this->nulldev_ << endl;
    for ( size_t i = 0; i < this->lambda.size(); i++){
        this->lambda[i] = 1.0 / (2.0+(double)i);
        LASSOgivenLambda currentLasso( x, wsaf, this->lambda[i], this->ju, this->g, this->ix);
        this->beta.push_back(currentLasso.beta);
        this->intercept.push_back(currentLasso.intercept);
        this->df.push_back(currentLasso.df);
        this->devRatio.push_back(currentLasso.devRatio);
    }
    cout<<"beta.size = "<< beta.size() <<endl;
}


void DEploidLASSO::productOfxy(vector < vector <double> > &x, vector < double > &wsaf){

    //10780 do 10781 j=1,ni
            //if(ju(j).eq.0)goto 10781
            //g(j)=abs(dot_product(y,x(:,j)))
    //10781 continue

    for (size_t j = 0; j < this->nVars_; j++){
        if (this->ju[j] != 0){
            double tmpSum = 0;
            for ( size_t i = 0; i < this->nObs_; i++ ){
                tmpSum += wsaf[i] * x[i][j];
            }
            this->g.push_back(std::abs(tmpSum));
        }
    }

    assert(this->g.size() == nVars);
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


void DEploidLASSO::initialization(size_t nLambda){
    this->lambda = vector <double> (nLambda, 0.0);
    this->intercept = vector <double> (nLambda, 0.0);
    this->devRatio = vector <double> (nLambda, 0.0);
    this->df = vector <int> (nLambda, 0);
}


LASSOgivenLambda::LASSOgivenLambda(vector < vector <double> > &x, vector < double > &wsaf, double lambda, vector < double > ju, vector <double> g, vector <double> &ix){
    //cout << lambda << endl;
    size_t nObs = x.size();
    // check nObs > 1
    size_t nVars = x[0].size();
    this->initialization(nObs, nVars);

    vector <double> vp(nVars, 1.0);

// ulam is user defined lambdas...
// no = number of observations
// ni = number of predictor variables

//c   parm = penalty member index (0 <= parm <= 1)
//c        = 0.0 => ridge
//c        = 1.0 => lasso
double beta = 1.0;
double bta=beta; // lasso part
double omb=1.0-bta; // ridge part
double rsq=0.0; //TODO, just replace as rsq0


    double alm0 = 0.0;
    double alm  = lambda;
    double flmin = 1.0; // this is defined in glmnet.R



                        //dem=alm*omb
                        //ab=alm*bta
                        //rsq0=rsq
                        //jz=1
                        //tlam=bta*(2.0*alm-alm0)

                  //10860 do 10861 k=1,ni
                              //if(ix(k).eq.1)goto 10861
                              //if(ju(k).eq.0)goto 10861
                              //if(g(k).gt.tlam*vp(k)) ix(k)=1
                        //10861 continue
                  //10862 continue

                double dem=alm*omb;
                  double      ab=alm*bta;
                    double    rsq0=rsq;
                     double   jz=1;
                     double   tlam=bta*(2.0*alm-alm0);

                  for ( size_t k = 0; k < nVars; k++){
                      if ((ix[k] != 1) & (ju[k] != 0)){
                         if ( g[k] > (tlam*vp[k]) ){
                             ix[k] = 1.0;
                         }
                      }
                  }



                  10870 continue
                  10871 continue
                        if(iz*jz.ne.0) go to 10360
                  10880 continue
                        nlp=nlp+1
                        dlx=0.0

                  10890 do 10891 k=1,ni
                              if(ix(k).eq.0)goto 10891
                              gk=dot_product(y,x(:,k))
                              ak=a(k)
                              u=gk+ak*xv(k)
                              v=abs(u)-vp(k)*ab
                              a(k)=0.0
                              if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*dem)))
                              if(a(k).eq.ak)goto 10891
                              if(mm(k) .ne. 0)goto 10911
                              nin=nin+1
                              if(nin.gt.nx)goto 10892
                              mm(k)=nin
                              ia(nin)=k
                        10911 continue
                              del=a(k)-ak
                              rsq=rsq+del*(2.0*gk-del*xv(k))
                              y=y-del*x(:,k)
                              dlx=max(xv(k)*del**2,dlx)
                        10891 continue
                  10892 continue

                        if(nin.gt.nx)goto 10872
                        if(dlx .ge. thr)goto 10931
                        ixx=0

                  10940 do 10941 k=1,ni
                        if(ix(k).eq.1)goto 10941
                        if(ju(k).eq.0)goto 10941
                        g(k)=abs(dot_product(y,x(:,k)))
                        if(g(k) .le. ab*vp(k))goto 10961
                        ix(k)=1
                        ixx=1
                  10961 continue
                  10941 continue
                  10942 continue

                        if(ixx.eq.1) go to 10880
                        goto 10872
                  10931 continue
                        if(nlp .le. maxit)goto 10981
                        jerr=-m
                        return
                  10981 continue
                  10360 continue
                        iz=1
                  10990 continue
                  10991 continue
                        nlp=nlp+1
                        dlx=0.0

                  11000 do 11001 l=1,nin
                              k=ia(l)
                              gk=dot_product(y,x(:,k))
                              ak=a(k)
                              u=gk+ak*xv(k)
                              v=abs(u)-vp(k)*ab
                              a(k)=0.0
                              if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*dem)))
                              if(a(k).eq.ak)goto 11001
                              del=a(k)-ak
                              rsq=rsq+del*(2.0*gk-del*xv(k))
                              y=y-del*x(:,k)
                              dlx=max(xv(k)*del**2,dlx)
                        11001 continue
                  11002 continue

                        if(dlx.lt.thr)goto 10992
                        if(nlp .le. maxit)goto 11021
                        jerr=-m
                        return
                  11021 continue
                              goto 10991
                  10992 continue
                        jz=0
                        goto 10871
                  10872 continue







}


void LASSOgivenLambda::initialization(size_t nObs, size_t nVars){
    this->maxIteration_ = 100000;
    this->thresh_ = 1e-7;
    this->absTol_ = 1e-4;
    this->nVars_ = nVars;
    this->dfmax_ = nVars_ + 1;
    //pmax = min(dfmax * 2+20, nvars)
    this->nObs_ = nObs;
    this->beta = vector <double> (nObs, 0.0);
    this->df = 0;
    this->devRatio = 0.0;
    this->intercept = 0.0;

    this->vp = vector <double> (nVars, 1.0);
    this->vq = vector <double> (nVars, 1.0); // For now ...


      //bta=beta
      //omb=1.0-bta
      //ix=0

      //if(flmin .ge. 1.0)goto 10771
            //eqs=max(eps,flmin)
            //alf=eqs**(1.0/(nlam-1))
//10771 continue

      //rsq=0.0
      //a=0.0
      //mm=0
      //nlp=0
      //nin=nlp
      //iz=0
      //mnl=min(mnlam,nlam)
      //alm=0.0



}
