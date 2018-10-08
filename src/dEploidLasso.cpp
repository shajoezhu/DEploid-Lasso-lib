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

#include "dEploidLasso.hpp"


DEploidLASSO::DEploidLASSO(const vector < vector <double> > &x,
    const vector < double > &wsaf, size_t nLambda) : Lasso(x, wsaf, nLambda) {
    this->setChoiceAt(nLambda-1);
    this->determineTheCutOff();
    this->shrinkThePanel();
}


void DEploidLASSO::determineTheCutOff() {
    if (this->lambda.size() == 0) {
        return;
    }

    for (size_t i = 0; i < (this->lambda.size()-1); i++) {
        double diff = (this->devRatio[i+1] - this->devRatio[i]) /
                      this->devRatio.back();
        if (diff < 0.001) {
            this->setChoiceAt(i);
            break;
        }
    }
}


void DEploidLASSO::shrinkThePanel() {
    if (this->lambda.size() == 0) {
        return;
    }

    vector <double> candidateBetas = this->beta[this->choiceAt()];
    for (size_t i = 0; i < candidateBetas.size(); i++) {
        if (candidateBetas[i] > 0) {
            choiceIdx.push_back(i);
            choiceBeta.push_back(candidateBetas[i]);
            // cout <<  choiceBeta.back() << endl;
        }
    }
}
