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
#include<string>
#include <fstream>      // std::ifstream

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
    cout<< "Matrix size = "<< x.size() << " " << x[0].size() << endl;
    cout<< "Vector length = " << wsaf.size() << endl;

    // Initialize
    this->nObs_ = wsaf.size();
    this->initialization();

    this->lambda = vector <double> (100, 0.0);
    for ( size_t i = 0; i < 100; i++){
        this->lambda[i] = 1.0 / (2.0+(double)i);
        // Given lambda[i], solve for currentBeta
        vector <double> currentBeta = this->solveBetaGivenLabmda( x, wsaf, this->lambda[i]);
        beta.push_back(currentBeta);
    }
    cout<<"beta.size = "<< beta.size() <<endl;
}


void DEploidLASSO::initialization(){
    this->maxIteration_ = 100;
    this->realTol_ = 1e-4;
    this->absTol_ = 1e-4;
}


vector < double > DEploidLASSO::solveBetaGivenLabmda(vector < vector <double> > &x, vector < double > &wsaf, double lambda){
    vector <double> ret;

    return ret;
}
