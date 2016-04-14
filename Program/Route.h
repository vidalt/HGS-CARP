/*  ---------------------------------------------------------------------- //
    Hybrid Genetic Search for Arc Routing Problems -- HGS-CARP
    Copyright (C) 2016 Thibaut VIDAL

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//  ---------------------------------------------------------------------- */

#ifndef ROUTE_H
#define ROUTE_H

#include "Params.h"
#include "SeqData.h"
#include "Noeud.h"
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <math.h>
class Noeud ;
using namespace std ;

class Route
{

private:

// Access to the parameters of the problem
Params * params ;

// Access to the associated individual
Individu * individu ;

public:

// Index of the route
int cour ;

// Day in which the route occurs
int day ;

// Depot node associated to the route
Noeud * depot ;

// Vehicle associated to the route
Vehicle * vehicle ;

// Is this route feasible or not
bool isFeasible ;

// The current cost of the route (updated after any route change and in the beginning of the LS).
double currentRouteCost ;

// Reverses the order of the route
void reverse () ;

// Update the data structures associated to this route
// If the flag "isForPrinting" is set to true, then the construction of the Seqdata 0_i will contain the information
// to remember the orientation of the visits (and not only the overall cost information)
void updateRouteData (bool isForPrinting) ;

// coutInsertionClient[i][p] stores the best insertion cost of client [i] with pattern [p] in this route
// The pattern information is due to the CARP specificity
vector < vector <double> > coutInsertionClient ;

// placeInsertionClient[i][p] stores the best place of insertion of client [i] with pattern [p] in this route
// The pattern information is due to the CARP specificity
vector < vector <Noeud *> > placeInsertionClient ;

// For each node, a bool saying if all moves involving this node and route have been tested without success
vector <bool> nodeAndRouteTested ;

// Reset the computation of all insertion values
void initiateInsertions();

// little debugging test
void testSeqDatas();

Route(void);

Route(int cour, Noeud * depot, Vehicle * vehicle, Params * params, Individu * indiv, int day);

~Route(void);
};

#endif
