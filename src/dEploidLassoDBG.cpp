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
#include <assert.h>     // assert
#include <limits>       // std::numeric_limits<double>::infinity();
#include <iomanip>      // std::setw




bool DEploidLASSO::print_normalized_struff(){
    dout<<"First 10 elements of normalized x"<<endl;
    for (size_t i = 0; i < nVars_; i++){
        vector<double> xx = this->standardized_x_transposed[i];
        for (size_t j = 0; j < xx.size(); j++){
            dout << std::setw(9) << xx[j] <<", ";
        }
    }
    dout<<endl;

    dout<<"First 10 elements of normalized y"<<endl;
    for (size_t i = 0; i < 10; i++){
        dout << std::setw(8) << this->standardized_y[i] <<", ";
    }
    dout<<endl;

    return(true);
}


bool DEploidLASSO::print_initial_gk(){
    dout << "Initial gk:"<<endl;
    for (double gk : g){
        dout << gk << ", ";
    }
    dout << endl;

    return(true);
}


bool DEploidLASSO::print_homogeneous_input(){
    // DEBUG MESSAGE
    dout << "Variables: ";
    for (size_t i = 0; i < this->nVars_; i++){
        if ( ju[i] == 0 ) {dout << i <<", ";}
    }
    dout << " are homogeneous vectors." << endl;

    return(true);
}
